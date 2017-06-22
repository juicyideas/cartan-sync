function testContradiction

n = 150;
d = 3;
% nze = 0;
nb = n/d;

m = rnd.cov(nb);
m = abs(m);
assert(all(eig(m)>0),'m is not PSD')

c_M = cell(nb,nb);
c_ev = cell(nb,nb);
for i=1:nb
  % Diagonal element
%   T = rnd.cov(d,1);
  mag = 10^(3*rand);
  T = mag*rand(d);
  [U,S,~] = svd(T);
  S = diag( diag(S) + m(i,i) ); % make sv_min = m_ii
  c_M{i,i} = U*S*U';
  
  for j=i+1:nb
    % Off-diagonal elements
%     T = rand(d);
    mag = 10^(3*rand);
    T = mag*rand(d);
    [U,~,V] = svd(T);
    S = diag( linspace(0,m(i,j),d) ); % make sv_max = m_ij
    c_M{i,j} = U*S*V';
    c_M{j,i} = c_M{i,j}';
  end
end
M = cell2mat(c_M);
M = symmetrize(M);

% Check if PSD
ev = eig(M)'

for i=1:nb
  for j=i:nb
    c_ev{i,j} = sort( eig(c_M{i,j}) );
    c_ev{j,i} = c_ev{i,j};
  end
end

% Compute m from M again
mSV = testBDD(M,d*ones(1,nb));
assert(norm(m-mSV,'fro')<1e-10,'m \neq mSV')

if any(ev<0)
  warning('Contradiction')
  keyboard
end

return
A = rnd.cov(n,2)


A = rnd.cov(n);
C = rnd.cov(n);
B = rand(n,n);
M = [A  B
     B' C];
mSV = testBDD(M,n*ones(1,2));

s = svd(mSV)';
assert(s(end)>=0)

end

function case2

% Generate random positive definite m
% m = rnd.cov(2);
m = [3 1; 1 0.5];
% m = [300 100; 100 0.5];
% m = [300 10; 10 0.4];
% All elements positive
m = abs(m);
assert(all(eig(m)>0),'m is not PSD')

% Define compliant block matrix M
% A
A = rnd.cov(d,1);
[U,S,V] = svd(A);
% s = diag(S);
% s(end) = m(1,1); % put m as smallest sing val
% S = diag(s);
S = diag( diag(S) + m(1,1) );
% S = diag( 100*diag(S) + m(1,1) );
A = U*S*U';
% B
B = rand(d);
[U,S,V] = svd(B);
s = diag(S);
% s(1) = m(1,2);
% S = diag(s);
% S = diag( s - max(s) + m(1,2) );
S = diag( linspace(0,m(1,2),d) );
B = U*S*V';
% C
C = rnd.cov(d,1);
[U,S,V] = svd(C);
% s = diag(S);
% s(end) = m(2,2); % put m as smallest sing val
% S = diag(s);
S = diag( 100*diag(S) + m(2,2) );
C = U*S*U';
% Build complete matrix
M = [A  B
     B' C];
M = symmetrize(M);

end