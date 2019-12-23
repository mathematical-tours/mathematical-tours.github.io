function [U,v,x2,y2] = FitRotation(x,y,x1,y1)

% minimized error
E = @(U,v)norm( U*[x y]'+v - [x1 y1]', 'fro');
% min_{U,v} sum_i |U*x_i+b - y_i|^2
% M = sum_i x_i y_i^T

m  = mean([x y])';
m1 = mean([x1 y1])';
%
M  = ([x y]'-m) * ([x1 y1]'-m1)';
U = inv(sqrtm(M'* M)) * M';
% U = inv(sqrtm(M* M')) * M;
% [G,~,H] = svd(M); U = G*H'; % SAME
v = m1 - U*m;
% 
Z = U*[x y]' + v;
x2 = Z(1,:)'; y2 = Z(2,:)';

end