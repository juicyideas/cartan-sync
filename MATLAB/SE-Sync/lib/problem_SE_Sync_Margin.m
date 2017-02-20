function problem = problem_SE_Sync_Margin(problem_data, d,n,p, use_Cholesky)
% function problem = problem_SE_Sync_Margin(problem_data, r)
%
% Returns a Manopt problem structure corresponding to the problem
%   ...
% [SO(3) synchronization with marginalized translational data]
%
% See also linearcost in Riemannian_staircase
%
% Nicolas Boumal, UCLouvain, March 4, 2014.

if nargin < 3
    use_Cholesky = true;
end

% We optimize over the manifold M := St(d, r)^N, the N-fold product of the
% (Stiefel) manifold of orthonormal d-frames in R^r.
manifold = stiefelstackedfactory(n, d, p);
problem.M = manifold;

% The cost
problem.cost = @(X) cost(X, problem_data, use_Cholesky);
  function [trQYtY, YQ] = cost(X, problem_data, use_Cholesky)
    % This function computes and returns the value of the objective function
    % 0.5*tr(Q Y^T Y).  Optionally, it returns the product YQ as the second
    % argument   
    Yt = X;
    YQ = Qproduct(Yt, problem_data, use_Cholesky)';
    trQYtY = 0.5*trace(YQ * Yt);
  end

% The Euclidean (classical) gradient w.r.t. Y
problem.egrad = @(X) egrad(X, problem_data, use_Cholesky);
  function G = egrad(X, problem_data, use_Cholesky)
    % This function computes and returns the value of the Euclidean gradient of
    % the objective function: nabla F(Y) = YQ.
    G = Qproduct(X, problem_data, use_Cholesky);
  end

% The Euclidean (classical) Hessian w.r.t. Y along the tangent vector Ydot
problem.ehess = @(X, Xdot) ehess(X, Xdot, problem_data, use_Cholesky);
  function H = ehess(X, Xdot, problem_data, use_Cholesky)
    % This function computes and returns the value of the Euclidean Hessian at
    % the point Y evaluated along the tangent direction Ydot.
    H = Qproduct(Xdot, problem_data, use_Cholesky);
  end

% keyboard
% Q = problem_data.ConLap + problem_data.Qt;
% problem.precon = @(X, H) manifold.proj(X, Q\H);
% problem.precon = @(X, H) preconditioner(X, H, problem_data);

do_Check = false;
if do_Check
  checkgradient( problem ); pause
  checkhessian(  problem ); pause
  keyboard
end

end