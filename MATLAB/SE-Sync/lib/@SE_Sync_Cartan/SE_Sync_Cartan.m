classdef SE_Sync_Cartan < SE_Sync
  
  methods
    
    function this = SE_Sync_Cartan(measurements, Manopt_opts, SE_Sync_opts)
      
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
          s_X0 = struct('R',xhat.R','t',xhat.t');
          ManiSEd = stiefelcartanstackedfactory(this.d,this.d,this.n);
          X0 = ManiSEd.toStacked(s_X0);
          % lift estimate to the correct dimension of Stiefel(p,d)
          X0 = this.liftEstimate(X0, this.opts.r0 - this.d);
          
          init_time = toc(init_time_start);
%           keyboard
          
        case {'chordal','chord'}
          fprintf('Computing chordal initialization...\n');
          init_time_start = tic();
          [R0,t0] = chordal_initialization(this.measurements);    
          % set as structure with column block-vectors
          % (a stiefelcartanstackedfactory point)
          % this estimate is in SE(d), so create temporary manifold struct
          s_X0 = struct('R',R0','t',t0');
          ManiSEd = stiefelcartanstackedfactory(this.d,this.d,this.n);
          X0 = ManiSEd.toStacked(s_X0);
          % lift estimate to the correct dimension of Stiefel(p,d)
          X0 = this.liftEstimate(X0, this.opts.r0 - this.d);
          
          init_time = toc(init_time_start);
          
        case {'origin-SO-Sync'}
          fprintf('Origin initialization followed by SO-Sync...\n');
          % create SO(d)-synchronization subproblem to initialize
          SO_Sync_opts = this.opts;
          SO_Sync_opts.init = 'origin';
          % TODO: Increase tolerance to finish earlier (and faster?)
          init_time_start = tic();
          subsync_problem = SO_Sync(this.measurements, this.Manopt_opts, SO_Sync_opts);
          % Run solver
          [SDPval, Xopt, xhat, Fxhat, SO_Sync_info] = run( subsync_problem );
          
          % tailor output for initialization of SE-Sync
          s_X0 = struct('R',xhat.R','t',xhat.t');
          ManiSEd = stiefelcartanstackedfactory(this.d,this.d,this.n);
          X0 = ManiSEd.toStacked(s_X0);
          % lift estimate to the correct dimension of Stiefel(p,d)
          X0 = this.liftEstimate(X0, this.opts.r0 - this.d);
          
          init_time = toc(init_time_start);
          
        case {'orig','origin'}
          fprintf('Using conventional origin in manifold as initialization...\n');
          init_time_start = tic();
          [R0,t0] = origin_initialization(this.measurements);    
          % set as structure with column block-vectors
          % (a stiefelcartanstackedfactory point)
          % this estimate is in SE(d), so create temporary manifold struct
          s_X0 = struct('R',R0','t',t0');
          ManiSEd = stiefelcartanstackedfactory(this.d,this.d,this.n);
          X0 = ManiSEd.toStacked(s_X0);
          % lift estimate to the correct dimension of Stiefel(p,d)
          X0 = this.liftEstimate(X0, this.opts.r0 - this.d);
          
          init_time = toc(init_time_start);
          
        case {'random','rand'}
          % Use randomly-sampled initialization
          fprintf('Randomly sampling an initial point on %s ...\n', this.Manopt_problem.M.name());
          init_time_start = tic();
          % Sample a random point on the manifold as an initial guess
          % TODO: Is Manopt set?
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
      
      % decide if using preconditioning from options
      usePrecon = this.opts.Precon;
      
      % call function
      problem = problem_SE_Sync_Cartan(this.problem_data, this.d,this.n,p, usePrecon);
      
      % Save Manopt structure
      this.Manopt_problem = problem;
      
%       % TODO: Move other place?
%       warning('We are changing Delta_bar for trustregions in Cartan')
%       this.Manopt_opts.Delta_bar = 10*problem.M.typicaldist();
      
    end
    
    % Other convenient tools
    function X_ = liftEstimate(this,X,dim_lift)
      X_ = horzcat(X, zeros((this.d+1) * this.n,dim_lift));
    end
    
    function c_Lambda = compute_Lambda(this, Xopt)
      %function Lambda_blocks = compute_Lambda(SE_Sync_Cartan, Xopt)
      %
      % Given an estimated local minimum Xopt for the (possibly lifted) relaxation,
      % this function computes and returns the block-diagonal elements of the
      % corresponding Lagrange multiplier
      
      % Copyright (C) 2017 by Jesus Briales
           
      % preallocate n dxd zero blocks
      % c_Lambda  - symmetric part of multipliers
      % c_ALambda - anti-symmetric part, this is the feasibility error
      [c_Lambda,c_A] = deal( repmat({zeros(this.d)},1,this.n) );
      
      % TODO: Work in progress...
      % compute MX product
      M = this.problem_data.ConLapT;
      manifold = this.Manopt_problem.M;
      MXopt = M*Xopt;
      
      for k = 1:this.n       
        B = manifold.get_Ri(MXopt,k)' * manifold.get_Ri(Xopt,k);
        Lam_k = .5 * (B + B');
        c_Lambda{k} = sparse(Lam_k); %force sparse for future matrices
        c_A{k} = B - Lam_k;
      end
      
      % extract translation component error
      err_t = manifold.get_t(MXopt);
      % NOTE: Error measures are not currently used since the minimum eig
      %       of the penalized matrix is a stronger test and it is
      %       necessary anyway for computing a descending direction
      % err_t should be near zero for optimal translation residues
      % c_A   should be near zero for optimal rotation residues
    end
    
    function [lambda_min, v] = min_eig_penalizedMat(this, c_Lambda, Xopt)
      
      % augment list of Lambda blocks with scalar zeros for t components
      % this construction provides Lambda blocks in the appropriate order
      % just calling c_Lambda{:}
      c_Lambda(2,:) = {0};;
      % build Lambda matrix
      Lambda = blkdiag( c_Lambda{:} );
      
      % call function implementing this
      [lambda_min, v] = Q_minus_Lambda_min_eig_Cartan(Lambda, this.problem_data, Xopt, this.opts.eig_comp_rel_tol);
         
    end
    
    function [xhat,Fxhat] = recover_solution( this, Xopt )
           
      % Round the solution
      disp('Rounding solution...');
      tic();
      [Xhat,singular_values] = round_solution_Cartan(Xopt, this.problem_data); % TRANSPOSED!
      solution_rounding_time = toc();
      fprintf('Elapsed computation time: %g seconds\n\n', solution_rounding_time);
      
      % Recover the optimal translational estimates
      disp('Recovering translational estimates... Not done here!');
      fprintf('Elapsed computation time: %g seconds\n\n', 0);


      if nargout > 1
        % Evaluate objective at the rounded solution
        Fxhat = this.Manopt_problem.cost( Xhat );
      end
      
      % Return solution in structure form
      % using Rosen's convention
      ManiSEd = stiefelcartanstackedfactory(this.d,this.d,this.n);
      xhat.R = ManiSEd.get_R(Xhat)';
      xhat.t = ManiSEd.get_t(Xhat)';
    end
    
    function info = analyze_Hessian( this, X )    
      
      % Get dimensionality of current estimate
      p = size(X,2);
      
%       % Setup problem
%       % Use the minimal vec/mat operators defined internally
%       prblm = this.Manopt_problem;
%       prblm.M.vec = prblm.M.minvec;
%       prblm.M.mat = prblm.M.minmat;
%       prblm.M.vecmatareisometries = @() true;
%       
%       % Check that unull is in the nullspace of Hess_f(X)[unull]
%       % get random vector in the nullspace
%       % there are p zero eigenvalues in the Hessian due to unull
%       unull = kron(ones(this.n,1),[zeros(this.d,1);1]) / sqrt(this.n); % normalized
%       Unull = unull*randn(1,p);
%       Q = this.problem_data.ConLapT;
%       fprintf('Unull is in Ker(Q): ||res|| = %g\n',norm(Q*Unull,'fro'))
%       hess_Unull = getHessian(prblm, X, Unull);
%       fprintf('Unull is in Ker(Hess_f(X)[U]): ||res|| = %g\n',norm(hess_Unull,'fro'))
%       % So Hess_f(X)[U] has ALWAYS at least p zero eigenvalues
      
      % Taylor our problem manifold to avoid any extra zeros
      % The manifold should actually be a quotient manifold
      % Anchor p extra dof in Euclidean components
      prblm = this.Manopt_problem;
      M = prblm.M;
      access = @(A,i) A(i);
      prblm.M.vec = @(X,u_mat) access( M.minvec(X,u_mat), 1:M.dim()-p );
      prblm.M.mat = @(X,u_vec) M.minmat(X,[u_vec;zeros(p,1)]);
      prblm.M.vecmatareisometries = @() true;
      prblm.M.dim = @() M.dim( ) - p; % remove anchored dof
      
      % Compute condition number with NO preconditioner
      [condNum,NZE,min_nz_eig,max_eig,min_eigs] = condest_Hessian(prblm, X, 'noprecon');
      info.noprec = setOutput(condNum,NZE,min_nz_eig,max_eig,min_eigs);
      
      % Now compute WITH preconditioner
      [condNum,NZE,min_nz_eig,max_eig,min_eigs] = condest_Hessian(prblm, X, 'precon');
      info.prec   = setOutput(condNum,NZE,min_nz_eig,max_eig,min_eigs);
%       keyboard
      
      function info = setOutput(condNum,NZE,min_nz_eig,max_eig,min_eigs)
        info.condNum = condNum;
        info.NZE = NZE;
        info.min_nz_eig = min_nz_eig;
        info.max_eig = max_eig;
        info.min_eigs = min_eigs;
      end
      
      
      
%       % Taylor our problem manifold to avoid any extra zeros
%       % The manifold should actually be a quotient manifold
%       % Anchor p extra dof in Euclidean components
%       access = @(A,i) A(i);
%       prblm.M.vec = @(X,u_mat) access( prblm.M.minvec(X,u_mat), 1:M.dim()-p );
%       prblm.M.mat = @(X,u_vec) prblm.M.minmat(X,[u_vec;zeros(p,1)]);
%       prblm.M.vecmatareisometries = @() true;
%       prblm.M.dim = @() prblm.M.dim( ) - p; % remove anchored dof
%       keyboard
% 
%       % setup anchored versions to remove nullspace of Connection Laplacian
%       % this removes p dof due to the last R^p Euclidean component
%       
%       [condNum,NZE,min_nz_eig,max_eig] = condest_Hessian(this.Manopt_problem, X);
      
    end
    
  end
end