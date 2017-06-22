function [mSV,dSV,mOff] = testBDD( M,dims )
% Check if the squared matrix M with input partition given by dims
% is Block Diagonally Dominant

s = size(M);
assert(s(1)==s(2),'Not square')
assert(issymmetric(M),'Not symmetric')
n = numel(dims);
% cA = mat2cell( M, d*ones(1,n),d*ones(1,n) );
offsets = cumsum([1 dims]);
cumdims = cumsum(dims);

% Compute and store required sub-matrix norms for BDD test
dSV = zeros(n,1);
% Smallest SV for diagonal blocks
for i=1:n
%   idxs_i = (i-1)*d+ (1:d);
  idxs_i = offsets(i) : cumdims(i);
%   Aii = cA{i,i};
  Aii = full( M(idxs_i,idxs_i) );
  minSV = inv( norm(inv(Aii),2) );
  dSV(i) = minSV;
end
% Largest SV for off-diagonal blocks
mOff = zeros(n,n);
for i=1:n-1
%   idxs_i = (i-1)*d+ (1:d);
  idxs_i = offsets(i) : cumdims(i);
  for j=i+1:n
%     idxs_j = (j-1)*d+ (1:d);
    idxs_j = offsets(j) : cumdims(j);
%     Aij = cA{i,j};
    Aij = full( M(idxs_i,idxs_j) );
    maxSV = norm(Aij,2);
    mOff(i,j) = maxSV;
    mOff(j,i) = maxSV;
  end
end
assert(issymmetric(mOff),'Not symmetric')

% Convert to sparse (do earlier for efficiency!?)
mOff = sparse(mOff);
mSV = diag(dSV) + mOff;

% Check diagonally dominance of usual scalar matrix
% Sum of row terms
rowSum = sum(mOff,2);
% Compare to diagonal
diff = dSV - rowSum;
% keyboard

if all(diff>=0)
  disp('Matrix is BDD')
else
  disp('Matrix is not BDD')
end

end