function [X, singular_values, determinants] = round_solution_Cartan(Xopt, problem_data)
%function [R, singular_values, determinants] = round_solution_Cartan(Yopt, problem_data)
%
% Given an element Yopt in St(d, r)^n, this function rounds Yopt to an
% element of SO(d)^n

% Copyright (C) 2016 by David M. Rosen

% Get problem dimensionality
% get dimensions of the rounded solution
d = problem_data.d;
n = problem_data.n;
p = size(Xopt, 2);

% The solution may be *contaminated* by components coming from
% the translation observability nullspace, so we remove those
% by projecting onto the range space (orthogonal to nullspace)
unull = kron(ones(n,1),[zeros(d,1);1]) / sqrt(n); % normalized
Xopt = Xopt - unull*(unull'*Xopt); % orthogonal projector

[U, Sigma, V] = svd(Xopt, 'econ');
singular_values = diag(Sigma)';

Sigma_d = Sigma(1:d, 1:d);  %Xi_d is the upper-left dxd submatrix of Xi
U_d = U(:,1:d);  %U_d contains the first d columns of U
X = U_d*Sigma_d;

% TODO: Assert this is a rank-d matrix, otherwise clever projections
%       should be used (see ICRA17)

determinants = zeros(1, problem_data.n);
ManiSEd = stiefelcartanstackedfactory(d,d,n);
for k = 1:problem_data.n
    determinants(k) = det( ManiSEd.get_Ri(X,k) );
end
ng0 = sum(determinants > 0);

reflector = diag([ones(1, problem_data.d - 1), -1]);  % Orthogonal matrix that we can use for reversing the orientations of the orthogonal matrix subblocks of R

if ng0 == 0
    % This solution converged to a reflection of the correct solution
    X = X*reflector;
    determinants = -determinants;
elseif ng0 < problem_data.n
    disp('WARNING: SOLUTION HAS INCONSISTENT ORIENTATIONS!');
    
    % If more than half of the determinants have negative sign, reverse
    % them
    if ng0 < problem_data.n / 2
        determinants = -determinants;
        X = X*reflector;
    end
end

% Finally, project each element of R to SO(d) using a retraction
% even if tight, for numerical stability
% TODO: This can be done more efficient for the {SE(d)}^n case
%       if using the specific manifold instead of the more generic
X = ManiSEd.retr(X, ManiSEd.zerovec());

end

