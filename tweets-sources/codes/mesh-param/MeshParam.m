%%
% Shows the evolution of a mesh parameterization.

addpath('../toolbox/');
rep = MkResRep();

boundname = 'circle';
boundname = 'square';

name = 'nefertiti';
[X,F] = read_off([name '.off']);
n = size(X,2);
clear options; options.name = name;


options.verb = 0;
B = compute_boundary(F, options);
p = length(B);

clf; hold on;
plot_mesh(X,F,options);
camlight;
shading 'faceted';
lighting flat;
plot3( X(1,B),X(2,B),X(3,B), 'r', 'LineWidth', 3 );
saveas(gcf, [rep name '-original.png']);

% Laplacian with cotan weights
[L,W,D,Di] = MeshLaplacian(X,F);

% fixed position for parameterization
p = length(B);
switch boundname
    case 'circle'
        t = pi*1.05 + linspace(0,2*pi(),p+1); t(p+1) = [];
        Z = [cos(t); sin(t)];
    case 'square'
        s = 5;
        B = B([s:end 1:s-1]);
        q = floor(p/4);
        t = (0:(q-1))/q;
        t1 = (0:(p-3*q-1))/(p-3*q);
        Z = [[t, t*0+1, 1-t, t1*0]; ...
            [t*0, t, t*0+1, 1-t1]];
        Z = 2*Z-1;
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
Y = (L1 \ R')'; 

clf; 
PlotMesh2d(F,Y);

% display the dynamic
% Y <- Proj( Y - tau*(Y - D^{-1} W Y) )
Y = randn(2,n)*.6;
niter = 50; 
tau = .9;
niter = 140; 
ndisp = floor(linspace(1,5,niter));
lw = 2; ms = 25; v = 1.02;
nvsg = 1;
for i=1:niter
    Y = Y' - tau*(Y' - Di*(W*Y') ); Y = Y';
    Y(:,B) = Z;
    if ndisp(i)==1 || mod(i,ndisp(i))==1
        clf;
        PlotMesh2d(F,Y);
        plot(Z(1,[1:end 1]), Z(2,[1:end 1]), 'r.-', 'LineWidth', lw, 'MarkerSize', ms);
        axis equal; axis([-v v -v v]); axis off;
        drawnow;
        saveas(gcf, [rep boundname '-' znum2str(nvsg,2) '.png']);
        nvsg = nvsg + 1;
    end
end


% AutoCrop(rep, [boundname '-']); 
