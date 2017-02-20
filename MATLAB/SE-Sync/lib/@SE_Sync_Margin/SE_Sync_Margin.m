classdef SE_Sync_Margin < SE_Sync
  
  methods
    
    function this = SE_Sync_Margin(measurements, Manopt_opts, SE_Sync_opts)
      
      % call parent constructor
      this@SE_Sync(measurements, Manopt_opts, SE_Sync_opts);
      
    end
    
    function construct_matrix_formulation(this)
      
      %% Construct problem data matrices from input
      fprintf('\n\nINITIALIZATION:\n\n');
      disp('Constructing auxiliary data matrices from raw measurements...');
      aux_time_start = tic();
      problem_data = construct_problem_data(this.measurements);
      this.info.mat_construct_times = toc(aux_time_start);
      fprintf('Auxiliary data matrix construction finished.  Elapsed computation time: %g seconds\n\n', this.info.mat_construct_times);
      % store problem data as property
      this.problem_data = problem_data;
    end
    
    
    function initialize( this, mode )
           
      %% INITIALIZATION
      
      switch mode
        case {'chordal-SO-Sync'}
          fprintf('Chordal initialization followed by SO-Sync...\n');
          % create SO(d)-synchronization subproblem to initialize
          SO_Sync_opts = this.opts;
          SO_Sync_opts.init = 'chordal';
          % TODO: Increase tolerance to finish earlier (and faster?)
          init_time_start = tic();
          subsync_problem = SO_Sync(this.measurements, this.Manopt_opts, SO_Sync_opts);
          % Run solver
          [SDPval, Xopt, xhat, Fxhat, SO_Sync_info] = run( subsync_problem );
          
          % tailor output for initialization of SE-Sync
          % Set as column block-vector (as stiefelstackedfactory)
          X0 = transpose( xhat.R );
          % Lift estimate to the correct dimension of Stiefel(p,d)
          X0 = liftEstimate(this,X0, this.opts.r0 - this.d);
          init_time = toc(init_time_start);
%           keyboard
        case {'chordal','chord'}
          fprintf('Computing chordal initialization...\n');
          init_time_start = tic();
          Rchordal = chordal_initialization(this.measurements);
          % Set as column block-vector (as stiefelstackedfactory)
          X0 = transpose( Rchordal );
          % Lift estimate to the correct dimension of Stiefel(p,d)
          X0 = liftEstimate(this,X0, this.opts.r0 - this.d);
          init_time = toc(init_time_start);
        case {'orig','origin'}
          fprintf('Using conventional origin in manifold as initialization...\n');
          init_time_start = tic();
          R0 = origin_initialization(this.measurements);
          % Set as column block-vector (as stiefelstackedfactory)
          X0 = transpose( R0 );
          % Lift estimate to the correct dimension of Stiefel(p,d)
          X0 = liftEstimate(this,X0, this.opts.r0 - this.d);
          init_time = toc(init_time_start);
        case {'random','rand'}
          % Use randomly-sampled initialization
%           fprintf('Randomly sampling an initial point on St(%d,%d)^%d ...\n', this.d, this.opts.r0, this.n);
          fprintf('Randomly sampling an initial point on %s ...\n', this.Manopt_problem.M.name());
          init_time_start = tic();
          % Sample a random point on the manifold as an initial guess
          X0 = this.Manopt_problem.M.rand();
          init_time = toc(init_time_start);
        otherwise
          error('Unknown initialization option')
      end
      fprintf('Elapsed init computation time: %g seconds\n', init_time);
      
      % Save initialization as property
      this.X0 = X0;
      
      this.info.init_time = init_time;
      
    end
    
    function setup_manopt_problem(this, p)
      % function problem = problem_marginSE(problem_data, r)
      %
      % Returns a Manopt problem structure corresponding to the problem
      %   ...
      % [SO(3) synchronization with marginalized translational data]
      %
      % See also linearcost in Riemannian_staircase
      %
      % Nicolas Boumal, UCLouvain, March 4, 2014.
      
      if nargin < 2
        p = this.opts.r0;
      end
      
      % Read options from problem configuration
      use_Cholesky = this.opts.Cholesky;
      
      % call function
      problem = problem_SE_Sync_Margin(this.problem_data, this.d,this.n,p, use_Cholesky);
      
      % Save Manopt structure
      this.Manopt_problem = problem;
      
    end
    
    % Other convenient tools
    function X_ = liftEstimate(this,X,dim_lift)
      X_ = horzcat(X, zeros(this.d * this.n,dim_lift));
    end
    
    function c_Lambda = compute_Lambda(this, Xopt)
      
      % call function in Rosen's library
      Lambda_blocks = compute_Lambda(Xopt', this.problem_data, this.opts.Cholesky);
      
      % set output in the common format for the abstracted SE-Sync
      c_Lambda = mat2cell( Lambda_blocks, this.d, this.d*ones(1,this.n) );
    end
    
    function [lambda_min, v] = min_eig_penalizedMat(this, c_Lambda, Xopt)
      
      % set data in the original format for Rosen's library
      Lambda = cell2mat(c_Lambda);
      Yopt = Xopt';
      
      % call function in Rosen's library
      [lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, this.problem_data, Yopt, this.opts.eig_comp_rel_tol, this.opts.Cholesky);
      
    end
    
    function [xhat,Fxhat] = recover_solution( this, Xopt )
      
      % Interface with Rosen's convention (row block-vector)
      Yopt = Xopt';
      
      % Round the solution
      disp('Rounding solution...');
      tic();
      Rhat = round_solution(Yopt, this.problem_data);
      solution_rounding_time = toc();
      fprintf('Elapsed computation time: %g seconds\n\n', solution_rounding_time);
      
      % Recover the optimal translational estimates
      disp('Recovering translational estimates...');
      tic();
      that = recover_translations(Rhat, this.problem_data);
      translation_recovery_time = toc();
      fprintf('Elapsed computation time: %g seconds\n\n', translation_recovery_time);
      
      % Return solution in structure form
      xhat.R = Rhat;
      xhat.t = that;
      
      if nargout > 1
        % Evaluate objective at the rounded solution
        Fxhat = evaluate_objective(Rhat, this.problem_data, this.opts.Cholesky);
      end
    end
    
    function info = analyze_Hessian( this, X )
      
      % Setup problem
      % Use the minimal vec/mat operators defined internally
      prblm = this.Manopt_problem;
      prblm.M.vec = prblm.M.minvec;
      prblm.M.mat = prblm.M.minmat;
      prblm.M.vecmatareisometries = @() true;
      
      [condNum,NZE,min_nz_eig,max_eig,min_eigs] = condest_Hessian(prblm, X, 'noprecon');
%       keyboard
      info.noprec = setOutput(condNum,NZE,min_nz_eig,max_eig,min_eigs);
      info.prec = []; % No preconditioning
      
      function info = setOutput(condNum,NZE,min_nz_eig,max_eig,min_eigs)
        info.condNum = condNum;
        info.NZE = NZE;
        info.min_nz_eig = min_nz_eig;
        info.max_eig = max_eig;
        info.min_eigs = min_eigs;
      end
      
    end
    
  end
  
  methods (Static)
    
    
  end
end