function [lambdas, V] = Q_minus_Lambda_min_eig_matrix(Lambda, Q, Xopt, tol, num_eigs)
%function [lambdas, V] = Q_minus_Lambda_min_eig_matrix(Lambda, problem_data, Xopt, tol, num_eigs)
%
% Given the Lagrange multiplier Lambda corresponding to a critical point
% Xopt of the low-rank Riemannian optimization problem, this function
% computes and returns the num_eigs smallest eigenvalues of the matrix 
% Q - Lambda, together with their corresponding eigenvectors V. Here 'tol' 
% refers to the relative tolerance of the minimum eigenvalue computation 
% using MATLAB's 'eigs' function.
% 
% Lambda here is the appropriate penalization matrix,
% from the Lagrangian relaxation (e.g. for SE(d)-sync or SO(d)-sync).

% Copyright (C) 2017 by Jesus Briales

if nargin < 4
    tol = 1e-5;  % default value
end

if nargin < 5
    num_eigs = 1;
end

% First, estimate the maximum eigenvalue of Q - Lambda
% (this should be its norm in the typical case, 
% as explained in http://math.stackexchange.com/q/603375/344323)
eigs_opts.issym = true;
eigs_opts.isreal = true;

% This function returns the product (Q - Lambda)*x
MminusLambda = Q - Lambda;

% eigs_opts.tol = 100*tol;  %This estimate needn't be particularly sharp...
% lambda_max = eigs(MminusLambda, 1, 'LA', eigs_opts)
% Faster? builtin Matlab function for 2-norm of sparse matrix (power iter)
lambda_max = normest( MminusLambda );


% We shift the spectrum of M - Lambda by adding lambda_max_est*I; this
% improves the condition number (to ~2 in the typical case)
% *without* perturbing any of the eigenspaces in Q - Lambda; this has the
% effect of producing MUCH faster convergence when running the Lanczos
% algorithm
%
MminusLambda_shifted = MminusLambda + lambda_max*speye(size(MminusLambda));

% Now compute the minimum eigenvalue of M - Lambda + lambda_max_est * I
eigs_opts.tol = tol;  %This should be a reasonably sharp estimate

if nargin >= 3
    % In the (typical) case that exactness holds, the minimum eigenvector
    % will be 0, with corresponding eigenvectors the columns of Xopt, so
    % we would like to use a guess close to this.  However, we would not
    % like to use /exactly/ this guess, because we know that (M - Lambda)X
    % = 0 implies that X lies in the null space of M - Lambda, and
    % therefore an iterative algorithm will get "stuck" if we start
    % /exactly/ there.  Therefore, we will "fuzz" this initial guess by
    % adding a randomly-sampled perturbation that is small in norm relative
    % to the first column of Xopt; this enables us to excite modes other
    % than Xopt itself (thereby escaping from this subspace in the 'bad
    % case' that Xopt is not the minimum eigenvalue subspace), while still
    % remaining close enough that we can converge to this answer quickly in
    % the 'good' case
    
%     relative_perturbation = .03;
    relative_perturbation = 3.0;
    % Why /sqrt(d) ?
%     eigs_opts.v0 = Xopt(:,1) + (relative_perturbation / sqrt(d))*randn(n*(d+1), 1);
    N = size(Xopt,1);
    eigs_opts.v0 = Xopt(:,1) + relative_perturbation*randn(N, 1);
end

[V, shifted_lambda_min] = eigs(MminusLambda_shifted, num_eigs, 'SA', eigs_opts);
lambdas = shifted_lambda_min - lambda_max*eye(num_eigs);

end

