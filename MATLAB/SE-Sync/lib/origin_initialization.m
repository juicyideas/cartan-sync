function [R0, t0] = origin_initialization(measurements)
% function [R0, t0] = origin_initialization(measurements)
% 
% OUTPUT: 
%  -Rchordal:  A d x dn block matrix Rchordal containing the origin
%     orientations.
%  -tchordal [optional]:  A d x n block matrix containing the origin
%     positions.

% Copyright (C) 2016 by David M. Rosen

d = length(measurements.t{1});
n = max(max(measurements.edges));

% rotations
R03 = repmat(eye(d),[1 1 n]);
R0 = reshape(R03,d,d*n);
% translations
t0 = zeros(d,n);

end

