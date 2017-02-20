function M = stiefelcartanstackedfactory(p,d,n)
% {Stiefel(p, d) x Re^p}^n,
% represented as a stacked factory:
%   transpose( [R1,t1,...,Rn,tn] )
% We use transpose to get a column block-vector, which is a more usual
% occurrence for the expressions in which this form appears.
% 
% Unless otherwise stated, when we write X we will refer to a point
% in the manifold such that X = transpose( [R1,t1,...,Rn,tn] ).
% This is a matrix submanifold in R^{(d+1)n x p}
% 
% Other useful versions of the point representation will be:
%   s_X = struct('R',[R1,...,Rn]','t',[t1,...,tn]')
%   X3D = cat(3,[R1,t1],...,[Rn,tn]) = cat(3,T1,...,Tn)
% Another useful notation will be the matrix Ti=[Ri,ti].
% 
% Note in all cases the familiarity with the special Euclidean group,
% as the manifold Stiefel(p, d) x Re^p can be seen as a lifted version
% of SE(d) to the p-th dimension arising in some relaxation problems.
% 
% For now, this class will be a wrapper to the simple productmanifold,
% but providing a stacked interface.
% In the future, we hope to implement this more efficiently

% Original author: Jesus Briales, Feb. 06, 2017.
% Contributors: 
% Change log:

if ~exist('n', 'var') || isempty(n)
  n = 1;
end

% create core product manifold and interface utilities
elements = struct();
elements.R = stiefelstackedfactory(n, d, p);
elements.t = euclideanfactory(     n,    p);
M_ = productmanifold(elements);

M.toStruct = @toStruct;
  function s_X = toStruct(X)
    % returns the point components as two fields:
    %  - R is a stiefelstackedfactory, dn x p
    %  - t is a euclideanfactory, n x p
    % this is the inverse function of toStacked
    s_X.R = X(all_R_idxs,:);
    s_X.t = X(all_t_idxs,:);
  end

M.toStacked = @toStacked;
  function X = toStacked(s_X)
    % returns the point components stacked in the default fashion
    % this is the inverse function of toStruct
    X = zeros(n*(d+1),p);
    X(all_R_idxs,:) = s_X.R;
    X(all_t_idxs,:) = s_X.t;
  end

% define standard Manopt fields
if n == 1
  M.name = @() sprintf('Stiefel cartan manifold St(%d, %d) x R^%d',p,d,p);
elseif n > 1
  M.name = @() sprintf('Power Stiefel cartan manifold {St(%d, %d) x R^%d}^%d, stacked', p,d,p,n);
else
  error('n must be an integer no less than 1.');
end

M.dim = M_.dim;

M.inner = @(x, d1, d2) d1(:).'*d2(:);

M.norm = @(x, d) norm(d(:));

M.dist = @(x, y) error('stiefelcartanstacked.dist not implemented yet.');

M.typicaldist = @() sqrt(M.dim());

M.proj = @(X, U) toStacked( M_.proj(toStruct(X),toStruct(U)) );
M.tangent = M.proj;

M.egrad2rgrad = M.proj;
M.ehess2rhess = @ehess2rhess;
  function rhess = ehess2rhess(X, egrad, ehess, Xdot)
    rhess = toStacked( M_.ehess2rhess( toStruct(X),...
                                       toStruct(egrad),...
                                       toStruct(ehess),...
                                       toStruct(Xdot) ) );
  end

M.retr = @retraction;
  function Y = retraction(X, U, t)
    if nargin < 3
      t = 1.0;
    end
    Y = toStacked( M_.retr(toStruct(X),toStruct(U),t) );
  end
M.exp  = @exponential;
  function Y = exponential(X, U, t)
    if nargin < 3
      t = 1.0;
    end
    Y = toStacked( M_.exp( toStruct(X),toStruct(U),t) );
  end

M.hash = @(X) ['z' hashmd5(X(:))];

M.rand = @() toStacked( M_.rand() );
M.randvec = @(X) toStacked( M_.randvec( toStruct(X) ) );

M.lincomb = @matrixlincomb;
    
M.zerovec = @(x) zeros(n*(d+1), p);
M.transp = @(x1, x2, d) M.proj(x2, d);
    
M.vec = @(x, u_mat) vec( transpose(u_mat) ); % vectorize transposed version
M.mat = @(x, u_vec) transpose( reshape(u_vec,[p,n*(d+1)]) );
M.vecmatareisometries = @() true;

M.minvec = @minvec;
M.minmat = @minmat;
  function u = minvec(X,U)
    
    s_X = toStruct(X);
    s_U = toStruct(U);
    % apply vectorization in Stiefel components
    u_R = elements.R.minvec(s_X.R,s_U.R);
    % apply (simple) vectorization in Euclidean components
%     u_t = elements.t.vec(transpose(s_X.t),s_U.t);
%     u_t = elements.t.vec(s_X.t,s_U.t);
    % Use OUR implementation (allows for flexible transposing)
    u_t = vec(transpose(s_U.t));
    % fuse together
    u = [u_R;u_t];
    
    % NOTE: When anchoring, note s_X.t is tall (vec does NOT return last
    % row as consecutive elements!)
% %     % NOTE2: Fixed for now transposing s_X.t
  end
  function U = minmat(X,u)
    
    s_X = toStruct(X);
    u_R = u(1:elements.R.dim());
    u_t = u(elements.R.dim()+1:end);
    % build mat in Stiefel components
    s_U.R = elements.R.minmat(s_X.R,u_R);
    % build (simple) mat in Euclidean components
%     s_U.t = transpose( elements.t.mat(s_X.t,u_t) );
%     s_U.t = elements.t.mat(s_X.t,u_t);
    % Use OUR implementation (allows for flexible transposing)
    s_U.t = reshape(u_t,p,n)';
    % fuse together
    U = toStacked(s_U);
  end


% create convenient tools for accessing manifold point data
% Define list of indeces for different subgroups of variables:
% all row-indeces in the stacked representation X\in R^{(d+1)n x p}
all_idxs = 1:(n*(d+1));
% array of all indeces for t components
all_t_idxs = (d+1) : (d+1) : (n*(d+1));
% cell list of groups of indeces for i-th t component
m_t_idxs   = all_t_idxs;
% array of all indeces for R components
all_R_idxs = setdiff(all_idxs,all_t_idxs);
% cell list of groups of indeces for i-th R component
m_R_idxs   = reshape(all_R_idxs,d,n);
m_T_idxs   = reshape(all_idxs,d+1,n);

% accessor functions
% Note: Stacked versions are column block-vectors (components are trans)
%       Component versions are straight (not transposed)
M.get_R  = @(X) X(all_R_idxs,:);
M.get_Ri = @(X,i) transpose( X(m_R_idxs(:,i),:) );
M.get_t  = @(X) X(all_t_idxs,:);
M.get_ti = @(X,i) transpose( X(m_t_idxs(:,i),:) );
M.get_Ti = @(X,i) transpose( X(m_T_idxs(:,i),:) );

end
