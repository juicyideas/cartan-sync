%function [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = SE_Sync(measurements, Manopt_opts, SE_Sync_opts, Y0)
%
% SE-Sync: A certifiably correct algorithm for synchronization over the
% special Euclidean group
%
%
% INPUTs:
%
% measurements:  A MATLAB struct containing the data describing the special
%   Euclidean synchronization problem (see eq. (11) in the paper for
%   details). Specifically, measurements must contain the following fields:
%   edges:  An (mx2)-dimensional matrix encoding the edges in the measurement
%     network; edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
%     that the states x_i are numbered sequentially as x_1, ... x_n.
%   R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
%   t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
%   kappa:  An m-dimensional cell array whose kth element gives the
%     precision of the rotational part of the kth measurement.
%   tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% Manopt_opts [optional]:  A MATLAB struct containing various options that
%       determine the behavior of Manopt's Riemannian truncated-Newton
%       trust-region method, which we use to solve instances of the
%       rank-restricted form of the semidefinite relaxation.  This struct
%       contains the following [optional] fields (among others, see the
%       Manopt documentation)
%   tolgradnorm:  Stopping criterion; norm tolerance for the Riemannian gradient
%   rel_func_tol:  An additional stopping criterion for the Manopt
%     solver.  Terminate whenever the relative decrease in function value
%     between subsequenct iterations is less than this value (in the range
%     (0,1) ).
%   maxinner:  Maximum number of Hessian-vector products to evaluate as part
%      of the truncated conjugate-gradient procedure used to compute update
%      steps.
%   miniter:  Minimum number of outer iterations (update steps).
%   maxiter:  Maximum number of outer iterations (update steps).
%   maxtime:  Maximum permissible elapsed computation time (in seconds).
%
% SE_Sync_opts [optional]:  A MATLAB struct determining the behavior of the
%       SE-Sync algorithm.  This struct contains the following [optional]
%       fields:
%   r0:  The initial value of the maximum-rank parameter r at which to
%      start the Riemannian Staircase
%   rmax:  The maximum value of the maximum-rank parameter r.
%   eig_comp_rel_tol:  Relative tolerance for the minimum-eigenvalue
%      computation needed to verify second-order optimality using MATLAB's
%      eigs command (typical values here are on the order of 10^-5)
%   min_eig_lower_bound:  Lower bound for the minimum eigenvalue in order to
%      consider the matrix Q - Lambda to be positive semidefinite.  Typical
%      values here should be small-magnitude negative numbers, e.g. -10^-4
%   Cholesky:  A Boolean value indicating whether to compute orthogonal
%      projections onto the cycle space of G using a cached Cholesky
%      factorization of Ared*Ared' or by applying an orthogonal (QR)
%      decomposition.  The former method may be faster on smaller problems,
%      but the latter is more numerically stable [default: false]
%   init:  A string specifying the initialization procedure to use if no
%      initial point Y0 is passed.  Options are 'chordal' or 'random'.  If
%      no option is specified, 'chordal' is used as a default
%
% Y0:  [Optional]  An initial point on the manifold St(d, r)^n at which to
%      initialize the first Riemannian optimization problem.  If this
%      parameter is not passed, a randomly-sampled point is used instead.
%
%
% OUTPUTS:
%
% SDPval:  The optimal value of the semidefinite relaxation
% Yopt:  A symmetric factor of an optimal solution Zopt = Yopt' * Yopt for
%      the semidefinite relaxation.
% xhat: A struct containing the estimate for the special Euclidean
%   synchronization problem.  It has the following two fields:
%   Rhat:  A d x dn matrix whose (dxd)-block elements give the rotational
%   state estimates.
%   that: a d x n matrix whose columsn give the translational state estimates.
% Fxhat:  The objective value of the rounded solution xhat.
%
% SE_Sync_info:  A MATLAB struct containing various possibly-interesting
%   bits of information about the execution of the SE-Sync algorithm.  The
%   fields are:
%   mat_contruct_times:  The elapsed time needed to construct the auxiliary
%     system matrices contained in 'problem_data'
%   init_time:  The elapsed time needed to compute the initial point for
%     the Riemannian Staircase.
%   optimization_times:  A vector containing the elapsed computation times
%     for solving the optimization problem at each level of the Riemannian
%     Staircase.
%   SDPLRvals:  A vector containing the optimal value of the optimization
%     problem solved at each level of the Riemannian Staircase
%   min_eig_times:  A vector containing the elapsed computation times for
%     performing the minimum-eigenvalue computation necessary to check for
%     optimality of Yopt as a solution of the SDP after solving each
%     Riemannian optimization problem to first-order.
%   min_eig_vals:  A vector containing the corresponding minimum
%      eigenvalues.
%   total_computation_time:  The elapsed computation time of the complete
%      SE-Sync algorithm
%   manopt_info:  The info struct returned by the Manopt solver for the
%      during its last execution (i.e. when solving the last explored level
%      the Riemannian Staircase).
%
% problem_data:  A MATLAB struct containing several auxiliary matrices
% constructed from the input measurements that are used internally
% throughout the SE-Sync algorithm.  Specifically, this struct contains the
% following fields:
%   n:  The number of group elements (poses) to estimate
%   m:  The number of relative measurements
%   d:  The dimension of the Euclidean space on which these group elements
%       act (generally d is 2 or 3).
%   LWtau:  The Laplacian for the translational weight graph W^tau.
%   ConLap:  The connection Laplacian for the set of rotational
%       measurements; see eq. (15) in the paper.
%   A:  An oriented incidence matrix for the directed graph of
%       measurements; see eq. (7) in the paper
%   Ared:  The reduced oriented incidence matrix obtained by removing the
%       final row of A.
%   L:  A sparse lower-triangular Cholesky factor of the reduced Laplacian
%       of the translational weight graph
%   T:  The sparse matrix of translational observations defined in eq. (24)
%       in the paper
%   Omega:  The diagonal matrix of translational measurement precisions;
%       defined in eq. (23).
%   V:  The sparse translational data matrix defined in eq. (16) in the
%       paper.

% Copyright (C) 2016 by David M. Rosen

classdef SE_Sync < handle
  % Handle to preserve applied changes and steps in the problem
  
  properties
    % Problem dimensions
    d,n,m
    % Data
    measurements
    problem_data    
    % Initial estimate
    X0
    % Options
    opts
    
    % Manopt interface
    Manopt_problem
    Manopt_opts
    
    % Information returned along execution
    info
  end
  
  methods
    
    function this = SE_Sync(measurements, Manopt_opts, SE_Sync_opts)
      %% SE_Sync Constructor
      % Setup data and options

      fprintf('\n\n========== SE-Sync ==========\n\n');
            
      % Store raw data:
      this.measurements = measurements;
      % Set dimensionality variables
      this.d = length(measurements.t{1});
      this.n = max(max(measurements.edges));
      this.m = size(measurements.edges, 1);
      
      %% INPUT PARSING
      
      % SE-Sync settings:
      fprintf('ALGORITHM SETTINGS:\n\n');
      
      if nargin < 3
        disp('Using default settings for SE-Sync:');
        SE_Sync_opts = struct;  % Create empty structure
      else
        disp('SE-Sync settings:');
      end
      
      if isfield(SE_Sync_opts, 'r0')
        fprintf(' Initial level of Riemannian Staircase: %d\n', SE_Sync_opts.r0);
      else
        SE_Sync_opts.r0 = 5;
        fprintf(' Setting initial level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.r0);
      end
      
      if isfield(SE_Sync_opts, 'rmax')
        fprintf(' Final level of Riemannian Staircase: %d\n', SE_Sync_opts.rmax);
      else
        SE_Sync_opts.rmax = 7;
        fprintf(' Setting final level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.rmax);
      end
      
      if isfield(SE_Sync_opts, 'eig_comp_rel_tol')
        fprintf(' Relative tolerance for minimum eigenvalue computation in test for positive semidefiniteness: %g\n', SE_Sync_opts.eig_comp_rel_tol);
      else
        SE_Sync_opts.eig_comp_rel_tol = 1e-5;
        fprintf(' Setting relative tolerance for minimum eigenvalue computation in test for positive semidefiniteness to: %g [default]\n', SE_Sync_opts.eig_comp_rel_tol);
      end
      
      if isfield(SE_Sync_opts, 'min_eig_lower_bound')
        fprintf(' Lower bound for minimum eigenvalue in test for positive semidefiniteness: %g\n', SE_Sync_opts.min_eig_lower_bound);
      else
        SE_Sync_opts.min_eig_lower_bound = -1e-3;
        fprintf(' Setting lower bound for minimum eigenvalue in test for positive semidefiniteness to: %g [default]\n', SE_Sync_opts.min_eig_lower_bound);
      end
      
      
      if ~isfield(SE_Sync_opts, 'Cholesky')
        fprintf(' Using QR decomposition to compute orthogonal projection [default]\n');
        SE_Sync_opts.Cholesky = false;
      else
        if SE_Sync_opts.Cholesky
          fprintf(' Using Cholesky decomposition to compute orthogonal projection\n');
        else
          fprintf(' Using QR decomposition to compute orthogonal projection\n');
        end
      end
      
      if ~isfield(SE_Sync_opts, 'init')
        fprintf(' Initialization method: chordal [default]\n');
        SE_Sync_opts.init = 'chordal';
      else
        switch SE_Sync_opts.init
          case {'chordal-SO-Sync'}
            fprintf(' Initialization method: chordal + SO-Sync\n');
          case {'chordal','chord'}
            fprintf(' Initialization method: chordal\n');
          case {'random','rand'}
            fprintf(' Initialization method: random\n');
          case {'origin-SO-Sync'}
            fprintf(' Initialization method: origin + SO-Sync\n');
          case {'origin','orig'}
            fprintf(' Initialization method: origin\n');
          otherwise
            error('Initialization option "%s" not recognized!  (Supported options are "chordal-SO-Sync", "chordal" or "random"\n', SE_Sync_opts.init);
        end
      end
      
      
      
      fprintf('\n');
      
      %% Manopt settings:
      
      if nargin < 2
        disp('Using default settings for Manopt:');
        Manopt_opts = struct;  % Create empty structure
      else
        disp('Manopt settings:');
      end
      
      if isfield(Manopt_opts, 'tolgradnorm')
        fprintf(' Stopping tolerance for norm of Riemannian gradient: %g\n', Manopt_opts.tolgradnorm);
      else
        Manopt_opts.tolgradnorm = 1e-2;
        fprintf(' Setting stopping tolerance for norm of Riemannian gradient to: %g [default]\n', Manopt_opts.tolgradnorm);
      end
      
      if isfield(Manopt_opts, 'rel_func_tol')
        fprintf(' Stopping tolerance for relative function decrease: %g\n', Manopt_opts.rel_func_tol);
      else
        Manopt_opts.rel_func_tol = 1e-5;
        fprintf(' Setting stopping tolerance for relative function decrease to: %g [default]\n', Manopt_opts.rel_func_tol);
      end
      
      if isfield(Manopt_opts, 'maxinner')
        fprintf(' Maximum number of Hessian-vector products to evaluate in each truncated Newton iteration: %d\n', Manopt_opts.maxinner);
      else
        Manopt_opts.maxinner = 500;
        fprintf(' Setting maximum number of Hessian-vector products to evaluate in each truncated Newton iteration to: %d [default]\n', Manopt_opts.maxinner);
      end
      
      if isfield(Manopt_opts, 'miniter')
        fprintf(' Minimum number of trust-region iterations: %d\n', Manopt_opts.miniter);
      else
        Manopt_opts.miniter = 1;
        fprintf(' Setting minimum number of trust-region iterations to: %d [default]\n', Manopt_opts.miniter);
      end
      
      if isfield(Manopt_opts, 'maxiter')
        fprintf(' Maximum number of trust-region iterations: %d\n', Manopt_opts.maxiter);
      else
        Manopt_opts.maxiter = 300;
        fprintf(' Setting maximum number of trust-region iterations to: %d [default]\n', Manopt_opts.maxiter);
      end
      
      if isfield(Manopt_opts, 'maxtime')
        fprintf(' Maximum permissible elapsed computation time [sec]: %g\n', Manopt_opts.maxtime);
      end
      
      
      % Check if a solver was explicitly supplied
      if(~isfield(Manopt_opts, 'solver'))
        % Use the trust-region solver by default
        Manopt_opts.solver = @trustregions;
      end
      solver_name = func2str(Manopt_opts.solver);
      if (~strcmp(solver_name, 'trustregions') && ~strcmp(solver_name, 'conjugategradient') && ~strcmp(solver_name, 'steepestdescent'))
        error(sprintf('Unrecognized Manopt solver: %s', solver_name));
      end
      fprintf('\nSolving Riemannian optimization problems using Manopt''s "%s" solver\n\n', solver_name);
      
      % Set additional stopping criterion for Manopt: stop if the relative
      % decrease in function value between successive iterates drops below the
      % threshold specified in SE_Sync_opts.relative_func_decrease_tol
      if(strcmp(solver_name, 'trustregions'))
        Manopt_opts.stopfun = @(manopt_problem, x, info, last) relative_func_decrease_stopfun(manopt_problem, x, info, last, Manopt_opts.rel_func_tol);
      end
      
      % Set additional statistics to be computed and stored along the alg
      Manopt_opts.statsfun = @mystatsfun;
      function stats = mystatsfun(problem, x, stats)
        % Compute the condition number of the Riemannian hessian
        % TODO
%         keyboard
%         lambdas = hessianspectrum(problem, x);
%         hist(lambdas) % margin case
%         stats.eigsHess = lambdas;
      end
      
      % Save options as properties
      this.opts = SE_Sync_opts;
      this.Manopt_opts  = Manopt_opts;
      
    end
    
    % Define headers for customizable blocks
    construct_matrix_formulation(this)
    setup_manopt_problem(this)
    initialize( this, mode )
    [SDPval, Xopt, xhat, Fxhat, SE_Sync_info] = solve(this) % COMMON BLOCK
    [xhat,Fxhat] = round_solution( this, Xopt )
    
    function [SDPval, Xopt, xhat, Fxhat, SE_Sync_info] = run( this )
      % batch execution of all blocks
      
      % Construct matrix formulation
      construct_matrix_formulation(this);
      % Setup structure for solving the problem with Manopt
      setup_manopt_problem(this);
      % Initialize the problem
      initialize( this, this.opts.init );
      % Run solver for SE-Sync
      [SDPval, Xopt, xhat, Fxhat, SE_Sync_info] = solve( this );
    end
        
    %% Headers to other convenient tools
    X_ = liftEstimate(X,dim_lift)
    % pad with extra zeros, equivalent to increasing dim in Stiefel mani
    
    Lambda_blocks = compute_Lambda(this, Xopt)
    % compute Lambdas (dual SDP solution) from nullspace condition
    % in a global critical point for tight cases
    
    [lambda_min, v] = min_eig_penalizedMat(this, c_Lambda, Xopt)
    % compute minimum eigenvalue (and corresponding eigenvector)
    % of the penalized matrix M-blkdiag(c_Lambda{:})
    % this may define if the critical point is a global minimum
    % exploiting the frequent tightness of the Lagrangian relaxation
  end

end