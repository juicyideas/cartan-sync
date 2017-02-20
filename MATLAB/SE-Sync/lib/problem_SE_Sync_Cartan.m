function problem = problem_SE_Sync_Cartan(problem_data, d,n,p, usePrecon)
% function setup_manopt_problem(this, p)
%
% Returns a Manopt problem structure corresponding to the problem
%   ...
% [SE(3) synchronization (without marginalization)]
%
% See also linearcost in Riemannian_staircase
%
% Nicolas Boumal, UCLouvain, March 4, 2014.

if nargin < 5
  usePrecon = true;
end

% We optimize over the manifold M := St(d, r)^N x {R^r}^N,
% product of the (Stiefel) manifold of orthonormal d-frames in R^r
% with augmented translation variables
manifold = stiefelcartanstackedfactory(p, d, n);
problem.M = manifold;

% % Test on spectrum of block-diagonal and block-off-diagonal components
% bdgmask = logical( kron(speye(n),true(d+1)) );
% offmask = ~bdgmask;
% M  = problem_data.ConLapT;
% Mb = problem_data.ConLapT;
% Mo = problem_data.ConLapT;
% Mb(offmask) = 0;
% Mo(bdgmask) = 0;
% eigs(M, 5,'LM')
% eigs(Mo,5,'LM')
% eigs(Mb,5,'LM')
% eigs(M, 5,'SM')
% eigs(Mo,5,'SM')
% eigs(Mb,5,'SM')
% sqMb = sqrtm(Mb);
% M_ = sqMb \ M / sqMb;
% M_ = Mb \ M / Mb;
% eigs(M_,5,'LM')
% eigs(M_,5,'SM')

% The cost
problem.cost = @(X) cost(X, problem_data);
  function [f, QX] = cost(X, data)
    % This function computes and returns the value of the objective function
    % 0.5*tr(Q X^T X).  Optionally, it returns the product XQ as the second
    % argument
    Q = data.ConLapT;
    QX = Q*X;
    f = 0.5*trace(X'*QX);
  end

% The Euclidean (classical) gradient w.r.t. X
problem.egrad = @(X) egrad(X, problem_data);
  function G = egrad(X, data)
    % This function computes and returns the value of the Euclidean gradient of
    % the objective function: nabla f(X) = QX.
    Q = data.ConLapT;
    df_dX = Q*X;
    G = df_dX;
  end

% The Euclidean (classical) Hessian w.r.t. X along the tangent vector Xdot
problem.ehess = @(X, Xdot) ehess(X, Xdot, problem_data);
  function H = ehess(X, Xdot, data)
    % This function computes and returns the value of the Euclidean Hessian at
    % the point X evaluated along the tangent direction Xdot.
    Q = data.ConLapT;
    ddf_U_dX = Q*Xdot;
    H = ddf_U_dX;
  end

% The preconditioner (approximately) solves the equation Hess_f(X)[U]=T
% problem.precon = @(X, T) preconditioner(X, T, this.problem_data);
% problem.precon = @(X, T) preconditioner_pcg(X, T, this.problem_data);
if usePrecon
  problem.precon = @(X, T) preconditioner_direct(X, T, problem_data);
end

% Test correctness of Manopt problem
do_Check = false;
if do_Check
  checkgradient( problem );
  keyboard
  checkhessian(  problem );
  keyboard
end

end

function U = preconditioner_direct( X,T, problem_data )

% profile clear
% profile on

% set manifold
d = problem_data.d;
n = problem_data.n;
p = size(X,2);
mani = stiefelcartanstackedfactory(p, d, n);

% set preconditioner data 
% NOTE: currently non-reduced matrix, see what happens
% NOTE: derivative without half factor (for now, as in Rosen)
Q = problem_data.ConLapT;
redQ = problem_data.redConLapT;
redL  = problem_data.Ls_redConLapT;
redLt = problem_data.Ls_redConLapTt;
% redS = problem_data.Ls_S;
s     = problem_data.Ls_s;
st(s) = 1:length(s); % https://es.mathworks.com/matlabcentral/answers/81070-how-to-find-the-reverse-of-a-permutation-vector#answer_90791
% redIL = problem_data.IL_redConLapT;
% redILt= problem_data.IL_redConLapTt;
% Compare number of nnz in complete and incomplete Chol
% nnz(redL), nnz(redIL)
% Tests for permutation vectors:
% I = speye(size(redL));
% S = I(s,s);
% St= I(st,st);
% Check correctness of Chol decomposition:
% norm( redM - (redS*redL)*(redL'*redS'), 'fro' )
% norm( redL*redL' - redM(s,s), 'fro' )
% keyboard

% project input vector to tangent space (unnecessary?)
T = mani.proj(X,T);

% compute (approximate) inverse for preconditioner
% use anchored matrix to avoid bad conditioning
redT = T(1:end-1,:);
% redPreconU = redM\redEtaU;
% Using incomplete Cholesky
% temp = redIL  \ redT;
% temp = redILt \ temp;
% redU = temp;
% Using complete Cholesky with pivoting
temp = redL  \ redT(s,:);
temp = redLt \ temp;
redU = temp(st,:);
U = [redU;zeros(1,p)];

% anchor mean translation component?
% NOTE: This seems to be essential for numerical stability of the precon
% U.t = U.t - ones(n,1)*mean(U.t,1);

% project preconditioned vector to tangent space again
U = mani.proj(X,U);

% profile viewer

% keyboard
end

%% Test: Use pcg as suggested by Bamdev
function U = preconditioner_pcg( X,T, problem_data )

profile clear
profile on

% set manifold
d = problem_data.d;
n = problem_data.n;
p = size(X.R,2);
mani = stiefel_euclidean_productfactory(p, d, n);

% set preconditioner data 
% NOTE: currently non-reduced matrix, see what happens
% NOTE: derivative without half factor (for now, as in Rosen)
M = problem_data.ConLapT;
redM = problem_data.redConLapT;
redL = problem_data.Ls_redConLapT;
redS = problem_data.Ls_S;
% check correctness
% norm( redM - (redS*redL)*(redL'*redS'), 'fro' )
% keyboard

% setup PCG
mytol = 1e-6;
% mymaxiter = 20;
mymaxiter = 500;
fun_pcg = @(x) compute_matrix_system(x,M,mani.catRtv(T),mani,X);
myprecon = [];
myprecon = @(eta) applyprecon(eta,redL,redS,mani,X); % any improvement?

% solve with PCG
[xsol,pcg_flag,RELRES,ITER,RESVEC] = pcg(fun_pcg, vec( mani.catRtv(T) ), mytol, mymaxiter, myprecon);
U = reshape(xsol,(d+1)*n,p);
U = mani.catRtv_inv( U );

profile viewer
% keyboard

% check the solution is actually tangent
fprintf('Residues before projT of solution U\n')
fprintf('Is tangent? %g\n', norm(mani.catRtv(mani.projN(X,U)),'fro') )
fprintf('Eq residue: %g\n', norm(mani.projmT(X,M*mani.catRtv(U)) - mani.catRtv(T),'fro') )
% project solution to tangent space (for numerical roundoff)
U = mani.projT(X,U);
fprintf('Residues after  projT of solution U\n')
fprintf('Is tangent? %g\n', norm(mani.catRtv(mani.projN(X,U)),'fro') )
fprintf('Eq residue: %g\n', norm(mani.projmT(X,M*mani.catRtv(U)) - mani.catRtv(T),'fro') )

% build function for preconjugated function: P^-1 * M * x
% boo = @(x) myprecon( fun_pcg(x) );
% analyze spectrum of the new preconditioned function
% eigs(boo,mani.nd_*mani.p,5,'LM')

keyboard
end

function lhsx = compute_matrix_system(x,M,T,mani,X)
U = reshape(x,mani.nd_,mani.p);
% project input vector to tangent space
U = mani.projmT(X,U);
% compute LHS
LHSx = M*U;
% project LHS to tangent space again
LHSx = mani.projmT(X,LHSx);
% vectorize for interface with PCG
lhsx = LHSx(:);

% disp(norm(lhsx - T(:)));
end
function Peta = applyprecon(eta,redL,redS,mani,X)
etaU = reshape(eta,mani.nd_,mani.p);
% project input vector to tangent space
etaU = mani.projmT(X,etaU);
% compute inverse for preconditioner
% use anchored matrix to avoid bad conditioning
redEtaU = etaU(1:end-1,:);
% redPreconU = redM\redEtaU;
redPreconU = redS * (redL' \ (redL \ (redS' * redEtaU)));
PreconU = [redPreconU;zeros(1,mani.p)];
% anchor mean translation component?
% NOTE: This seems to be essential for numerical stability of the precon
% U.t = U.t - ones(n,1)*mean(U.t,1);

% project preconditioned vector to tangent space again
PreconU = mani.projmT(X,PreconU);
% vectorize for interface with PCG
Peta = PreconU(:);
end

%% Test: Solving preconditioner equation (matrix form)
function U = preconditioner( X,T, problem_data )
% U = preconditioner( X, T )
% Solve the equation Hess_f(X)[U]=T (inverse of the Hessian),
% where both U and T lie in T_X_Mani (tangent space).
% We use the simplified Hessian Proj_X(d^2f_dX^2) as an approximation.

profile clear
profile on

% set manifold
d = problem_data.d;
n = problem_data.n;
p = size(X.R,2);
mani = stiefel_euclidean_productfactory(p, d, n);

% set preconditioner data 
% NOTE: currently non-reduced matrix, see what happens
% NOTE: derivative without half factor (for now)
M = problem_data.ConLapT;
% extract current point in the manifold (as blocks)
[c_R,c_t] = mani.toCell(X);

% set all matrices appearing in the formulation
% NOTE: Make matrices sparse
shuffle_dd = sparse( build_vecProj(d,d) );
Pp = speye(d^2) + shuffle_dd; % projection matrix
inv_Pp = sparse( round( 4*pinv(full(Pp)) ) / 4 ); % quick hack, but exact
Pm = speye(d^2) - shuffle_dd; % projection matrix
I_kr_Pp = kron(speye(n),Pp);
I_kr_invPp = kron(speye(n),inv_Pp);
I_kr_Pm = kron(speye(n),Pm);
% block-diagonal matrix with current point
c_I_kr_R = cell(1,n);
for i=1:n
  c_I_kr_R{i} = [ kron(speye(d), c_R{i}); zeros(p,d^2) ];
end
B_X = blkdiag(c_I_kr_R{:});
% inverse of data matrix M inside kronecker product
% using Cholesky decomposition M = S*Lp*Lp'*S'
Lp = problem_data.Ls_redConLapT; % Lower triangular factor
S  = problem_data.Ls_S; % Permutation matrix in Chol decomp
% L = problem_data.Ls_S * problem_data.Ls_redConLapT; % TODO: IMPORTANT keep permutation out, to have triangular matrix!!!!!!!!!!
L = problem_data.L_redConLapT; % TODO: IMPORTANT keep permutation out, to have triangular matrix!!!!!!!!!!
% L(:,end+1) = 0;
% L(end+1,:) = 0;
krL = kron(L,speye(p)); % Using direct chol
krS = kron(S, speye(p));
krLp = kron(Lp,speye(p));  % Using chol with pivoting (requires perm with S too)
krLt = krLp';
redB_X = B_X;
redB_X(end-p+1:end,:) = [];
% redB_X(dropIdxs,:) = [];
% Recompute symmetric factors with Cholesky
% redLB_X = krL\redB_X;
% build projection matrix from vec(S_i) to svec(S_i)
Ps   = sparse(build_vec2vecsProj(d));
% c_Ps = repmat( {Ps}, 1,n );
% I_kr_Ps = blkdiag(c_Ps{:});
I_kr_Ps = kron(speye(n),Ps); % USE THIS

% solve linear system
useChol = true;
mT = mani.catRtv(T);
red_mT = mT(1:end-1,:);
if useChol
%   A = +I_kr_Ps * I_kr_Pp * (redLB_X' * redLB_X) * I_kr_Ps';
%   A = +I_kr_Ps * I_kr_Pp * redLB_X' * redLB_X * I_kr_Ps';
%   Al = +I_kr_Ps * I_kr_Pp * redLB_X';
%   Ar = krL \ redB_X * I_kr_Ps';
%   A = Al*Ar;
  invLvecT = krLp\ (krS'*vec(red_mT'));
  b = -I_kr_Ps * I_kr_Pp * redB_X' * krS * ( krLp' \ invLvecT );
else
  % Use complete matrices, padded with zeros
%   A = +I_kr_Ps * I_kr_Pp * B_X' * inv_krM * B_X * I_kr_Ps';
%   b = -I_kr_Ps * I_kr_Pp * B_X' * inv_krM * vT;
  % Use reduced matrices, in concrete inv_krM uses inv_redM
  A = +I_kr_Ps * I_kr_Pp * redB_X' * inv_krredM * redB_X * I_kr_Ps';
  b = -I_kr_Ps * I_kr_Pp * redB_X' * inv_krredM * vec(red_mT');
  % Debug:
%   norm( redB_X' * inv_krredM * redB_X - B_X' * inv_krM * B_X, 'fro' )
end
% s = A \ b;
tic
% y = Al \ b;
y = krLp' * krS' * redB_X * I_kr_invPp * I_kr_Ps' * b;

% Make chain of operations here too

fprintf('Al \\ b (underdetermined) time: %g\n',toc)
% tic
% s = Ar \ y;
% fprintf('Ar \\ y (underdetermined) time: %g\n',toc)
% tic
% [Q,R,E] = qr(Ar,0);
% fprintf('Aqr(Ar,0)                time: %g\n',toc)
% foo = Ar'*Ar;
% svds( foo, 5, 'largest' )
% svds( foo, 5, 'smallest' )
% condNum = svds(foo,1,'largest') / svds(foo,1,'smallestnz');
% fprintf('Original condition number: %f\n',condNum)
% Try to create simple preconditioner
% Precon = I_kr_Ps * redB_X' * krL * krL' * redB_X * I_kr_Ps';
% Precon = I_kr_Ps * redB_X' * krL * krL' * redB_X * I_kr_Ps';
krM = kron(problem_data.redConLapT,speye(p));
Precon = I_kr_Ps * redB_X' * krM * redB_X * I_kr_Ps';
Precon = symmetrize(Precon);
% condNum = svds(foo*Precon,1,'largest') / svds(foo*Precon,1,'smallestnz')
% fprintf('New condition number: %f\n',condNum)

% % % Check L*L' is the same as original M
% % foo = L*L';
% % foo( abs(foo)<1e-6 ) = 0;
% % spy( foo )
% % % Check if inv of foo as some pattern (as M)
% % foo = Ar'*Ar;
% % foo( abs(foo)<1e-6 ) = 0;
% % spy( foo )
% % inv_foo = inv( foo );
% % norm( foo * inv_foo - eye(750), 'fro' )
% % norm( inv_foo - Precon, 'fro' )
% % inv_foo( abs(inv_foo)<1e-3 ) = 0;
% % norm( foo * inv_foo - eye(750), 'fro' )
% % norm( inv_foo - Precon, 'fro' )
% % inv_foo( abs(inv_foo)<1e-0 ) = 0;
% % norm( foo * inv_foo - eye(750), 'fro' )
% % norm( inv_foo - Precon, 'fro' )
% % inv_foo( abs(inv_foo)<10 ) = 0;
% % norm( foo * inv_foo - eye(750), 'fro' )
% % norm( foo * Precon - eye(750), 'fro' )
% % norm( inv_foo - Precon, 'fro' )
% % subplot(1,2,1), spy(Precon), subplot(1,2,2), spy( inv_foo )
% % %
% % surf( inv_foo ), shading flat, colorbar
% % coo = inv_foo * foo;
% % coo( abs(coo)<1e-6 ) = 0;
% % spy( coo )
% % spy( inv_foo )
% % spy( Precon )
% % %
% % boo = inv( L'\inv(L) );
% % boo( abs(boo)<1e-6 ) = 0;
% % spy( problem_data.redConLapT )
% % spy( boo )
% % norm( problem_data.redConLapT - boo, 'fro' )
% % %
% % s = (Ar'*Ar) \ Ar'*y; % normal equation



% compare number of nz elements: it's preferable to not to compute Ar'*Ar
% nnz(redB_X)+nnz(I_kr_Ps)+nnz(krL), nnz(redB_X)+nnz(I_kr_Ps)+nnz(krL) nnz(Ar), nnz(foo)

  function Mx = fun_Mx(x)
    % Mx = m(x)
    % Return result of the product Ar'*Ar*x
    temp = I_kr_Ps' * x;
    temp = redB_X  * temp;
    temp = krLp \ (krS'*temp);
    temp = krS * (krLt \ temp);
    temp = redB_X' * temp;
    Mx   = I_kr_Ps * temp;
  end
% Test function
% r = randn(size(foo,1),1);
% norm( foo*r - fun_Mx(r) )
tic
% b_y = Ar'*y;
% b_y = krL \ redB_X * I_kr_Ps' * y;
b_y = I_kr_Ps * (redB_X' * krS * (krLp' \ y));
[s,pcg_flag,RELRES,ITER,RESVEC] = pcg(@fun_Mx,b_y);
% [s,pcg_flag,RELRES,ITER,RESVEC] = pcg(@fun_Mx,b_y, [],[], Precon);
% [s,pcg_flag,RELRES,ITER,RESVEC] = pcg(@fun_Mx,b_y, [],50, Precon);
% [s,pcg_flag,RELRES,ITER,RESVEC] = pcg(@fun_Mx,b_y, 1e-3, 10);
% [s,pcg_flag,RELRES,ITER,RESVEC] = pcg(@fun_Mx,b_y, 1e-5, 1000);
% [s,pcg_flag,RELRES,ITER,RESVEC] = pcg(@fun_Mx,b_y, 1e-5, 1000, Precon);
fprintf('pcg method #iters %d, time: %g\n',ITER,toc)
% Check step by step
% norm( Ar \ y - s )
% norm( Ar*s-y )
% norm( redB_X * I_kr_Ps'*s - krL*y )
% y_ = redB_X' * krL * y;
% s = I_kr_Ps' \ y_;
% s = I_kr_Ps * redB_X' * krL * y;
% keyboard
% s = Ar \ (Al \ b);
% rebuild symmetric matrices
c_S = mat2cell( reshape( I_kr_Ps'*s, d,d*n ), d, d*ones(1,n) );
% solve normal vectors
c_N = cell(1,n);
for i=1:n
  c_N{i} = [c_R{i}*c_S{i}, zeros(p,1)];
end
N = cell2mat(c_N)';
N = mani.catRtv_inv(N);
% solve tangent vector
U = M \ (mani.catRtv(N) + mani.catRtv(T));
U = mani.catRtv_inv(U);
% anchor mean translation component
% NOTE: This seems to be essential for numerical stability of the precon
U.t = U.t - ones(n,1)*mean(U.t,1);
% [c_UR,c_Ut] = mani.toCell(U);

% Project result to have true tangent vector as output
Uproj = mani.proj(X,U);
Nproj = M * mani.catRtv(Uproj) - mani.catRtv(T);
Nproj = mani.catRtv_inv(Nproj);

%% check solution
doCheck = true;
if doCheck
  % check if all the conditions in the equation are solved
  fprintf('Preconditioner tests:\n')
  % tangent condition error
  tangentError = mani.catRtv(mani.proj(X,U)) - mani.catRtv(U);
  fprintf('Tangent error: %E\n',norm(tangentError,'fro'))
  % normal error
  normalError  = mani.catRtv(mani.proj(X,N));
  fprintf('Normal  error: %E\n',norm(normalError,'fro'))
  % equation error
  MU = mani.catRtv_inv( M*mani.catRtv(U) );
  eqError = mani.catRtv(mani.proj(X,MU)) - mani.catRtv(T);
  fprintf('Eq error: %E\n',norm(eqError,'fro'))
  % Note: if numerical error is big it can be due to bad conditioning of M
  
  
  % check if all the conditions in the equation are solved
  fprintf('Preconditioner tests (after projection):\n')
  % tangent condition error
  tangentError = mani.catRtv(mani.proj(X,Uproj)) - mani.catRtv(Uproj);
  fprintf('Tangent error: %E\n',norm(tangentError,'fro'))
  % normal error
  normalError  = mani.catRtv(mani.proj(X,Nproj));
  fprintf('Normal  error: %E\n',norm(normalError,'fro'))
  % equation error
  MUproj = mani.catRtv_inv( M*mani.catRtv(Uproj) );
  eqError = mani.catRtv(mani.proj(X,MUproj)) - mani.catRtv(T);
  fprintf('Eq error: %E\n',norm(eqError,'fro'))
  % Note: if numerical error is big it can be due to bad conditioning of M
  
  
  profile viewer
  keyboard
end

% Return true tangent vector
U = Uproj;

end