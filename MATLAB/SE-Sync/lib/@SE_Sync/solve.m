function [SDPval, Xopt, xhat, Fxhat, SE_Sync_info] = solve(this)

% The maximum number of levels in the Riemannian Staircase that we will
% need to explore
max_num_iters = this.opts.rmax - this.opts.r0 + 1;

% Allocate storage for state traces
optimization_times = zeros(1, max_num_iters);
SDPLRvals = zeros(1, max_num_iters);
min_eig_times = zeros(1, max_num_iters);
min_eig_vals = zeros(1, max_num_iters);

% Counter to keep track of how many iterations of the Riemannian Staircase
% have been performed
iter = 0;

% Set initialization to that previously stored in the problem
X0 = this.X0;

%%  RIEMANNIAN STAIRCASE
for r = this.opts.r0 : this.opts.rmax
  iter = iter + 1;  % Increment iteration number
  
  % Starting at X0, use Manopt's truncated-Newton trust-region method
  % to descend to a first-order critical point.
  
  fprintf('RIEMANNIAN STAIRCASE (level r = %d):\n', r);
%   profile clear
%   profile on
  [Xopt, Fval, manopt_info, this.Manopt_opts] = manoptsolve(this.Manopt_problem, X0, this.Manopt_opts);  % complete
%   profile viewer
%   keyboard
  SDPLRval = Fval(end);
  
  % Store the optimal value and the elapsed computation time
  SDPLRvals(iter) = SDPLRval;
  optimization_times(iter) = manopt_info(end).time;
  
  % Augment Xopt to next Stiefel dimension by padding with zeros; this
  % preserves Xopt's first-order criticality while ensuring that it is
  % rank-deficient
  Xplus = this.liftEstimate(Xopt,1);
  
  fprintf('\nChecking second-order optimality...\n');
  % At this point, Xplus is a rank-deficient critical point, so check
  % 2nd-order optimality conditions
  
  % Compute Lagrange multiplier matrix Lambda corresponding to Xplus
  c_Lambda = this.compute_Lambda(Xopt);
  
  % Compute minimum eigenvalue/eigenvector pair for Q - Lambda
  tic();
  [lambda_min, v] = this.min_eig_penalizedMat(c_Lambda, Xopt);
  min_eig_comp_time = toc();
  
  % Store the minimum eigenvalue and elapsed computation times
  min_eig_vals(iter) = lambda_min;
  min_eig_times(iter) = min_eig_comp_time;
  
  if( lambda_min > this.opts.min_eig_lower_bound)
    % Yopt is a second-order critical point
    fprintf('Found second-order critical point! (minimum eigenvalue = %g, elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time);
    break;
  else
    fprintf('Saddle point detected (minimum eigenvalue = %g,  elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time);
    % lambda_min is a negative eigenvalue of Q - Lambda, so the KKT
    % conditions for the semidefinite relaxation are not satisfied;
    % this implies that Xplus is a saddle point of the rank-restricted
    % semidefinite optimization.  Fortunately, the eigenvector v
    % corresponding to lambda_min can be used to provide a descent
    % direction from this saddle point, as described in Theorem 3.9 of
    % the paper "A Riemannian Low-Rank Method for Optimization over
    % Semidefinite Matrices with Block-Diagonal Constraints".
    
    % Define the vectorXYdot := v * e_{r+1}; this is tangent to the
    % lifted manifold at Xplus and provides a direction of
    % negative curvature
    disp('Computing escape direction...');
    Xdot = horzcat( this.Manopt_problem.M.zerovec(), v);
    
    % Compute the directional derivative of F at Xplus along Xdot
    dF0 = trace( Xdot'*this.Manopt_problem.egrad(Xplus) );
    if dF0 > 0
      Xdot = -Xdot;
    end
    
    % Augment the dimensionality of the Stiefel manifolds in
    % preparation for the next iteration
    setup_manopt_problem(this, r+1)
    
    % Perform line search along the escape direction Xdot to escape the
    % saddle point and obtain the initial iterate for the next level in
    % the Staircase
    disp('Line searching along escape direction to escape saddle point...');
    tic();
    [stepsize, X0] = linesearch_decrease(this.Manopt_problem, Xplus, Xdot, SDPLRval);
    line_search_time = toc();
    fprintf('Line search completed (elapsed computation time %g seconds)\n', line_search_time);
  end
end

fprintf('\n\n===== END RIEMANNIAN STAIRCASE =====\n\n');

% Return optimal value of the SDP (in the case that a rank-deficient,
% second-order critical point is obtained, this is equal to the optimum
% value obtained from the Riemannian optimization
SDPval = SDPLRval;

% POST-PROCESSING
% Recover solution in the original domain
[xhat,Fxhat] = this.recover_solution( Xopt );

fprintf('Value of SDP solution F(Y): %g\n', SDPval);
fprintf('Norm of Riemannian gradient grad F(Y): %g\n', manopt_info(end).gradnorm);
fprintf('Value of rounded pose estimate xhat: %g\n', Fxhat);
fprintf('Suboptimality bound of recovered estimate: %g\n\n', Fxhat - SDPval);

%       fprintf('Total elapsed computation time: %g seconds\n\n', total_computation_time);

% Output info
SE_Sync_info.SDPLRvals = SDPLRvals(1:iter);
SE_Sync_info.mat_construct_times = this.info.mat_construct_times;
SE_Sync_info.init_time = this.info.init_time;
SE_Sync_info.optimization_times = optimization_times(1:iter);
SE_Sync_info.min_eig_vals = min_eig_vals(1:iter);
SE_Sync_info.min_eig_times = min_eig_times(1:iter);
SE_Sync_info.manopt = manopt_info;
SE_Sync_info.total_computation_time = []; %TODO: Sum as dep prop

fprintf('\n===== END SE-SYNC =====\n');

end