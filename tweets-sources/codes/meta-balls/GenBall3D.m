function [F,R,T] = GenBall3D(c,w,n, phi, Tp)

% c: center
% w: weight
% n: size of image

K = length(w);

col = distinguishable_colors(K+1); 
if K>=4
    col(4,:) = [];
end
mm = @(x)mod(x+1/2,1)-1/2; % in [-1/2,1/2] modulo 1
mm = @(x)x;


t = linspace(0,1,n);
[X,Y,Z] = ndgrid(t,t,t);

resh = @(u)repmat(reshape(u, [1 1 1 K]), [n n n]);
D = mm( repmat(X,[1 1 1 K]) - resh(c(:,1)) ).^2 + ...
    mm( repmat(Y,[1 1 1 K]) - resh(c(:,2)) ).^2 + ...
    mm( repmat(Z,[1 1 1 K]) - resh(c(:,3)) ).^2;
f = phi(D);
F = sum(f .* resh(w), 4);


fC = f .* resh(max(w,0));
FC = sum(fC, 4);

% colors
W = fC ./ repmat(FC, [1 1 1 K]);
R = zeros(n,n,n,3);
for i=1:3
    for k=1:K
        R(:,:,:,i) = R(:,:,:,i) + W(:,:,:,k) * col(k,i); 
    end
end
% masking
a = sort(abs(F(:)), 'descend');
T = a(round(abs(Tp)*end));
% M = F>T;
% RF = R .* repmat(M, [1 1 1 3]) + repmat(1-M, [1 1 1 3]);

end