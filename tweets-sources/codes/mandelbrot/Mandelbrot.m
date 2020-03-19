%%
% Plot Mandelbrot as complex polynomials

addpath('../toolbox/');
rep = MkResRep();


xmax = 3; 
ymax = 3;
n = 512;
y = linspace(-1.6,1.6,n);
x = linspace(-2.2,1.2,n);
[Y,X] = meshgrid(y,x);

Z = X + 1i*Y;
P = Z;
pmax = 8;
smax = 8;
it = 0; 
for p=1:pmax
    P1 = P.^2 + Z;
	m = 1; 
    Q = Z.^(2^(p-1));  Q1 = Z.^(2^(p));
    for s=1:smax
        t = (s-1)/smax;
        clf;
        PolyDisp((1-t.^m)*P./Q+t.^m*P1./Q1,x,y);
        % PolyDisp( P.^(t+1) + t*Z ,x,y);
        % PolyDisp(  + t*Z ,x,y);
        axis equal; axis square;
        drawnow;
        it = it+1;
        %saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);
    end
    P = P.^2 + Z;
end