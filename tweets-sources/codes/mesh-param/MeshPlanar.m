%%
% Test for parameterization of planar mesh on a disk


name = 'cat';
name = 'elephant';
name = 'france';

%%
% Shows the evolution of a mesh parameterization.

addpath('../toolbox/');
addpath('../mesh-2d/');
rep = MkResRep(name);
libpath();

%% load from a binary image

f = load_image(name, 512);
f = sum(f,3); f = double((f/max(f(:)))>.5);
[C,h] = contour(f,[.5,.5]);
m = C(2,1);  c = C(:,2:m+1); c = c(1,:)'+1i*c(2,:)';
c = (c-1)/(max(size(f))-1);
node0 = [real(c), imag(c)];

% curve re-sampling
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );

% re-sample
p = 400;
c = node0(:,1) + 1i*node0(:,2);
c = resample(c,p);
node = [real(c), imag(c)];
n = size(node,1);
edge = [1:n; [2:n 1]]';


% generate
opt.disp = Inf;
opt.rho2 = 1.00;
% opt.size1 = 1.33;
[X,etri, F,tnum] = refine2(node,edge,[],opt); % hfun
% polish
[X,etri,F,tnum] = smooth2(X,etri,F,tnum,opt);
X = X'; F = F';
n = size(X,2);
X = [X;zeros(1,n)];
X = X - repmat(mean(X,2), [1 n]);
X = X/mean(sqrt(sum(X.^2,1)));



options.verb = 0;
B = compute_boundary(F, options);
p = length(B);


VC = ones(n,1);
clf; hold on;
h = patch('vertices',X','faces',F','FaceVertexCData',VC, 'FaceColor', 'interp');
plot3( X(1,B),X(2,B),X(3,B), 'r', 'LineWidth', 3 );
% shading interp;
axis off;


% Laplacian with cotan weights
[L,W,D,Di] = MeshLaplacian(X,F);

% fixed position for parameterization
d = Inf;
Z = 0;
for rho = linspace(0,2*pi,100)
    p = length(B);
    t = linspace(0,2*pi(),p+1)+rho; t(p+1) = []; % pi*1.05 +
    z = [cos(t); sin(t)];
    if norm(X(1:2,B)-z, 'fro')<d
        Z = z;
        d = norm(X(1:2,B)-z, 'fro');
    end
    if norm(X(1:2,B)-z(:,end:-1:1), 'fro')<d
        Z = z(:,end:-1:1);
        d = norm(X(1:2,B)-z(:,end:-1:1), 'fro');
    end
end


% fixed positions Laplacian
L1 = L;
L1(B,:) = 0;
for i=1:length(B)
    L1(B(i),B(i)) = 1;
end
% RHS
R = zeros(2, n);
R(:,B) = Z;

% Solve linear system
Y0 = (L1 \ R')'; 
Y = [Y0;zeros(1,n)];


% Display as texture
r = 800; % texture precision
x = linspace(-1,1,r);
[V,U] = meshgrid(x,x);
a = cos( 12*pi* sqrt(U.^2+V.^2) );
b = cos( 10.5*pi* atan2(V,U) );
I = double(  a.*b>0 );
VT = (Y0+1)/2;
opts.EdgeColor = 'none';
opts.PSize = 128;

q = 50; 
for i=1:q
    t = (i-1)/(q-1);
    Xt = (1-t)*X + t*Y;
    clf;
    patcht(F',Xt',F',VT',I, opts);
    axis off; axis ij; hold on;
    plot3( Xt(1,B([1:end 1])),Xt(2,B([1:end 1])),Xt(3,B([1:end 1])), 'r', 'LineWidth', 3 );
    axis off; axis ij;
    axis equal; axis tight;
    % axis([-1 1 -1 1]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, 'anim-');

