%%
% Displays the bassin of attraction of Newton method on the complex plane.

if not(exist('cnt'))
    cnt = 1;
end

addpath('../toolbox/');
rep = MkResRep(['interp/' num2str(cnt) '/']);

if not(exist('theta'))
    % rotated newton
    theta = 1.4; rho = .7;
    % basic newton
    rho = 0;
end

% click to generate the roots of the polynomials

P = {};
for k=1:2
    P{k} = [];
    clf; hold on;
    while true
        axis equal; axis([0 1 0 1]);
        if k==2
            plot(real(P{1}), imag(P{1}), 'b.', 'MarkerSize', 25);
        end
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if button==3
            break;
        end
        P{k}(end+1) = a+1i*b;
    end
    Pc{k} = poly(P{k});
end

% pad with zero
d = max(length(Pc{1}), length(Pc{2}));
for k=1:2
    Pc{k} = [zeros(d-length(Pc{k}),1); Pc{k}(:)];
end

% grid size
N = 512;
% width of the window
A = 1/2;
% center
u = 1/2+1i/2; 

nimg = 50;  % #images
Pt_Old = P{1}; Pt = P{1}; 
for i = 1:nimg
    t = (i-1)/(nimg-1);  
    % interpolate polynomial
    Pct = (1-t)*Pc{1}+t*Pc{2};
    % roots
    PtNew = roots(Pct);
    % re-align the roots with the previous one.
    PtNew = PtNew(randperm(length(PtNew)));
    [Pt,Pt_Old] = deal(RealignRoots(PtNew,Pt_Old),Pt);    
    % evaluation function and derivative 
    f = @(Z)polyval(Pct, Z);
    tau = 1e-6;
    Df = @(Z)(f(Z+tau)-f(Z))/tau;
    % Newton map
    Phi = @(z)z - ( 1 - rho*exp(1i*theta) ) * f(z)./Df(z);       
    % generate image
    [C1,x,y] = PlotNewtonFractal(Phi,Pt,N,A,u);
    % imwrite(C1, [rep 'newton-fractal-' num2str(q) '.png'], 'png');
    % imwrite(C1, [rep 'newton-' znum2str(i,2) '.png'], 'png');
    %
    clf; hold on;
    imagesc(x,y,permute(C1, [2 1 3]));
    plot(real(Pt), imag(Pt), 'k.',  'MarkerSize', 25, 'Color', 'k');
    axis image; axis([0 1 0 1]); axis off;
    drawnow;
    saveas(gcf, [rep 'newton-' znum2str(i,2) '.png']);
end

cnt = cnt+1;

% AutoCrop(rep, 'newton-');