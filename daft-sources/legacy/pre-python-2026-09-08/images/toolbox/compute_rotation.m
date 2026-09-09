function M = compute_rotation(v,alpha)

% compute the 3D rotation matrix around vector v
%   and with angle alpha.
%
%   M = compute_rotation(v,alpha)
%
%   Copyright (c) 2004 Gabriel Peyré

v = reshape(v,3,1);
v = v/norm(v,'fro');

S = [0   -v(3) v(2);
    v(3)   0  -v(1);
    -v(2)  v(1) 0];
M = v*transp(v) + cos(alpha)*(eye(3) - v*transp(v)) + sin(alpha)*S;
