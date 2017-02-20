function problem_data = construct_problem_data(measurements)
%function problem_data = construct_problem_data(measurements)
%
% This helper function accepts a MATLAB struct containing the raw 
% measurements specifying a special Euclidean synchronization problem, and 
% returns another struct containing the data matrices required by the 
% SE-Sync algorithm.  Formally:
%
% INPUT: A MATLAB struct 'measurements' containing the following fields
% (see eq. (11) in the long-form version of the paper for details):
% edges:  An (mx2)-dimension encoding the edges in the measurement network;
%     edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform from pose i to pose j.  NB:  This indexing scheme 
%     requires that the states x_i are numbered sequentially as 
%     x_1, ... x_n.
% R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
% t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
% kappa:  An m-dimensional cell array whose kth element gives the precision
%     of the rotational part of the kth measurement. 
% tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% 
%
% OUTPUT:  A MATLAB struct 'problem_data' containing various data matrices
% that the SE-Sync algorithm requires:
%
% d:  The dimensional parameter of the special Euclidean group over which 
%     the estimation takes place (typically d = 2 or 3).
% n:  The number of states to be estimated.
% m:  The number of available measurements.
% LWtau:  The Laplacian for the translational weight graph W^tau.
% ConLap:  The connection Laplacian for the rotational measurements (see
%     eq. (15) in the paper.
% A:  The oriented incidence matrix for the underlying directed graph of
%     measurements (see eq. (7) in the paper).
% Ared:  The reduced incidence matrix obtained by removing the final 
%     row of A.
% L:  The lower-triangular Cholesky factor of the reduced translational
%     weight graph Laplacian.
% Omega:  The diagonal matrix of translational matrix precisions (see eq.
%     (23) in the paper).
% T:  The sparse matrix of translational observations definedin equation 
%     (24) in the paper.
% V:  The matrix of translational observations defined in equation 
%     (16) in the paper

%Given the 'measurements' struct returned by 'readG2oDataset3D', this
%function constructs and returns the data matrices defining the pose-graph
%relaxation problem

% Copyright (C) 2016 by David M. Rosen

% Set additional variables
d = length(measurements.t{1});
n = max(max(measurements.edges));
m = size(measurements.edges, 1);
problem_data = struct('d',d,'n',n,'m',m);


% Construct the oriented incidence matrix for the underlying directed
% connection graph of SE(d) measurements
% Each vertex corresponds to a block [Ri,ti]
tic();
[AT,OmegaT] = construct_connection_incidence_matrix(measurements);
problem_data.ConInc   = AT;
problem_data.ConOmega = OmegaT;
t = toc();
fprintf('Constructed connection incidence and weight matrices in %g seconds\n', t);

% Construct the connection Laplacian for the SE(d) measurements
% Construct also its anchored version (drop last column and row)
tic();
problem_data.ConLapT = symmetrize( AT * OmegaT * AT' ); % symmetric by construction, fix numerical issues
problem_data.redConLapT = problem_data.ConLapT(1:end-1,1:end-1);
t = toc();
fprintf('Constructed connection Laplacian in SE(d) in %g seconds\n', t);

% Construct the Cholesky factor for the reduced SE Connection Laplacian
% tic();
% problem_data.L_redConLapT = chol(problem_data.redConLapT, 'lower');
% t = toc();
% fprintf('Computed lower-triangular factor of reduced connection Laplacian in SE(d) in %g seconds\n', t);
% tic();

% Construct the Cholesky factor for the reduced SE Connection Laplacian
% using pivoting to produce sparser results
% [Ls,p,S] = chol(problem_data.redConLapT, 'lower','matrix');
[Lp,cholFlag,s] = chol(problem_data.redConLapT, 'lower','vector');
assert(cholFlag==0,'Error computing Cholesky')
problem_data.Ls_redConLapT  = Lp;
problem_data.Ls_redConLapTt = Lp';
% problem_data.Ls_S = S;
problem_data.Ls_s = s;
% norm( S*Ls*Ls'*S'-problem_data.redConLapT, 'fro' )
% subplot(1,2,1), spy(problem_data.redConLapT)
% subplot(1,2,2), spy(S*Ls*Ls'*S')
t = toc();
fprintf('Computed lower-triangular factor with pivoting of reduced connection Laplacian in SE(d) in %g seconds\n', t);
% keyboard

% Set indeces to access sections of interest in the matrix
% Define list of indeces for different subgroups of variables:
% all row-indeces in the stacked representation X\in R^{(d+1)n x p}
% There are two different sets of indeces involved: Vertex and Edge kind
all_vIdxs   = 1:(n*(d+1));
all_eIdxs   = 1:(m*(d+1));
% array of all indeces for t components
all_t_vIdxs = (d+1) : (d+1) : (n*(d+1));
all_t_eIdxs = (d+1) : (d+1) : (m*(d+1));
% cell list of groups of indeces for i-th t component
% m_t_idxs   = all_t_idxs;
% array of all indeces for R components
all_R_vIdxs = setdiff(all_vIdxs,all_t_vIdxs);
all_R_eIdxs = setdiff(all_eIdxs,all_t_eIdxs);
% cell list of groups of indeces for i-th R component
% m_R_idxs   = reshape(all_R_idxs,d,n);
% m_T_idxs   = reshape(all_idxs,d+1,n);


% Construct connection Laplacian for the rotational measurements
tic();
% problem_data.ConLap = construct_connection_Laplacian(measurements);
% Take only part due to rotation observations
AR = AT(all_R_vIdxs,all_R_eIdxs);
OmegaR = OmegaT(all_R_eIdxs,all_R_eIdxs);
problem_data.ConLap = symmetrize( AR * OmegaR * AR' ); % symmetric by construction, fix numerical issues
% problem_data.ConLap = problem_data.ConLapT(all_R_idxs,all_R_idxs); % THIS IS WRONG!
t = toc();
fprintf('Constructed rotational connection Laplacian in %g seconds\n', t);

% Construct the Cholesky factor for the SO-Connection Laplacian
% using pivoting to produce sparser results
% Lp stands for Lower matrix with Pivoting
[Lp,cholFlag,s] = chol(problem_data.ConLap, 'lower','vector');
assert(cholFlag==0,'Error computing Cholesky of SO-Conn Laplacian')
problem_data.Lp_ConLapR   = Lp;
problem_data.Lp_ConLapRt  = Lp';
problem_data.Lp_ConLapR_s = s;
% norm( S*Ls*Ls'*S'-problem_data.redConLapT, 'fro' )
% subplot(1,2,1), spy(problem_data.redConLapT)
% subplot(1,2,2), spy(S*Ls*Ls'*S')
t = toc();
fprintf('Computed lower-triangular factor with pivoting of SO(d)-Connection Laplacian in %g seconds\n', t);
% keyboard

% Construct the oriented incidence matrix for the underlying directed graph
% of measurements
tic();
problem_data.A = construct_incidence_matrix(measurements);
t = toc();
fprintf('Constructed oriented incidence matrix in %g seconds\n', t);

% Construct the reduced oriented incidence matrix
problem_data.Ared = problem_data.A(1:problem_data.n-1, :);

tic();
[T, Omega] = construct_translational_matrices(measurements);
V = construct_V_matrix(measurements);
t = toc();
fprintf('Constructed translational observation and measurement precision matrices in %g seconds\n', t);

problem_data.T = T;
problem_data.Omega = Omega;
problem_data.V = V;


% Construct the Laplacian for the translational weight graph
tic();
LWtau = problem_data.A * problem_data.Omega * problem_data.A';
t = toc();
fprintf('Constructed Laplacian for the translational weight graph in %g seconds\n', t);
problem_data.LWtau = LWtau;


% Construct the Cholesky factor for the reduced translational weight graph
% Laplacian
tic();
problem_data.L  = chol(LWtau(1:end-1, 1:end-1), 'lower');
problem_data.Lt = transpose( problem_data.L ); % store transpose for speed
t = toc();
fprintf('Computed lower-triangular factor of reduced translational weight graph Laplacian in %g seconds\n', t);

% Cache a couple of various useful products
fprintf('Caching additional product matrices ... \n');

problem_data.sqrt_Omega = spdiags(sqrt(spdiags(Omega)), 0, problem_data.m, problem_data.m);
problem_data.sqrt_Omega_AredT = problem_data.sqrt_Omega * problem_data.Ared';
problem_data.sqrt_Omega_T = problem_data.sqrt_Omega * problem_data.T;

% Compute Incomplete Cholesky Factorization of Connection Laplacian
% tic();
% opts = struct();
% opts.type = 'nofill';
% opts.michol = 'on';
% opts.type = 'ict';
% opts.droptol = 1e-3;
% IL = ichol(problem_data.redConLapT,opts);
% IL = ichol(problem_data.ConLapT,opts);
% t = toc();
% fprintf('Computed Incomplete Cholesky of reduced connection Laplacian in SE(d) in %g seconds\n', t);
% problem_data.IL_redConLapT  = IL;
% problem_data.IL_redConLapTt = IL';
% keyboard

end