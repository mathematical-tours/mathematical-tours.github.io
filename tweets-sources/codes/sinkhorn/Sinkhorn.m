%%
% Sinkhorn iterations in the simplest bistochastic case.

name = 'transcale';
name = 'folding';
addpath('../toolbox/');
rep = MkResRep(name);

n = 32; 


% piecewise constant upsampling
S = @(k,p)floor(1:1/k:p+1-1/k);
Upsc = @(x,k)x(S(k,size(x,1)),S(k,size(x,2)));
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% 
x = linspace(-1,1,n);
switch name
    case 'translation'
        s = .1; d = .8;
        y = .8 + x;
    case 'folding'
        s = .3;
        y = 2*abs(x)-3;
        y = 5*abs(x);
    case 'transcale'
        s = .3;
        y = 1+4*x;
end

[Y,X] = meshgrid(y,x);
K = exp( -(X-Y).^2/(2*s^2) );


K = K/sum(K(:));

niter = 50; 
a = ones(n,1); b = ones(n,1);
for i=1:niter
    %
    a = 1./(K*b);
    m1 = b.*(K'*a);
    A1 = diag(a)*K*diag(b);
    %
    b = 1./(K'*a);
    m2 = a.*(K*b);
    A2 = diag(a)*K*diag(b);
    %    
    clf; 
    subplot(2,2,1);
    imagesc(-A1); axis image; axis off;
    subplot(2,2,2);
    imagesc(-A2); axis image; axis off;
    subplot(2,2,3);
    plot(m1); 
    axis([1 n 0 3]);
    subplot(2,2,4); hold on;
    plot([1 n], [1 1], 'k--', 'LineWidth', 1);
    plot(m2, 'LineWidth', 2);
    axis([1 n 0 2]); box on;
    drawnow;
    %
    clf; hold on;
    plot([1 n], [1 1], 'k--', 'LineWidth', 1);
    plot(m2, 'LineWidth', 2);
    axis([1 n 0 2]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    SetAR(1/2);
    saveas(gcf, [rep 'marginal-' znum2str(i,2) '.png'], 'png');
    % rendering
    u = 4; 
    AA = Upsc(rescale(-A1),u);
    t = (i-1)/(niter-1);
    c = [t 0 1-t];
    B = zeros(n*u,n*u,3);
    for s=1:3
        B(:,:,s) = (1-AA)*c(s) + AA*1;
    end
    imwrite(B, [rep 'coupling-' znum2str(i,2) '.png'], 'png');    
end


% AutoCrop(rep, 'marginal-');