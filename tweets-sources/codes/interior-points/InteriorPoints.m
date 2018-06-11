%%
% Display interior point for linear programming

if not(exist('cnt'))
    cnt = 1;
end

addpath('../toolbox/');
rep = MkResRep(num2str(cnt));

% Generate a polygon
% if not(exist('Z'))
it = 0;
clf; hold on;
Z = [];
while true
    axis([0 1 0 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(end+1) = a+1i*b;
end
% end
n = length(Z);

% Generate inequalities
b = zeros(n,1);
A = zeros(n,2);
for k=1:n
    u = Z(k); v = Z(mod(k,n)+1);
    m = exp(1i*pi/2)*( u-v ); m = m/abs(m);
    A(k,:)=[real(m), imag(m)];
    % <m,a>=b
    b(k) = real(u*conj(m));
end
b = b(:);

% display polygon
p = 256*2;
%
t = linspace(0,1,p);
[Y,X] = meshgrid(t,t);
V = [X(:)'; Y(:)'];

% polygon A*x<b
U = repmat(b(:),[1 p*p]) - A*V;
Q = prod(U>=0);
D = min(U,[],1);

resh = @(Q)reshape(Q,p,p);
clf; imagesc(resh(Q));

% function to minimize 
r = randn(2,1);
f0 = resh(V'*r);

% pre-compute path
q = 50;  % #images
epslist = 10.^(-linspace(3,-1,q));
% pre-compute path
W = [];
for i=1:q
    epsilon = epslist(i);
    f = resh( r'*V - epsilon*sum( log(max(U,1e-40)) ) ); f(Q==0) = Inf;
    [~,m] = min(f(:)); [x,y] = ind2sub([p p], m); 
    W(i) = t(x) + 1i * t(y);
end

% render
for i=1:q
    epsilon = epslist(i);
    f = resh( r'*V - epsilon*sum( log(max(U,1e-20)) ) );
    %
    J = find(D>.01);
    vmin = min(f(J));
    vmax = max(f(J));
    %
    clf;
    RenderPolygon(Q,f,vmin,vmax,Z);
    plot(W, 'r', 'LineWidth', 2);
    plot(W(i), 'r.', 'MarkerSize', 25);
    drawnow;
    saveas(gcf, [rep 'ip-' znum2str(i,2) '.png']);
end
cnt = cnt+1;
% AutoCrop(rep, 'ip-');

    