%%
% Display noise property.

% Remove low freq bar
% Keep global data and make global windowing

name = 'poisson-disk';
name = 'poisson';
name = 'square';
name = 'hexagonal';

addpath('../toolbox/');
rep = MkResRep(name);

% size for discrete fourier
n = 1024*2;


sp = .008*n; % spacing for poisson disk

% target number
X = poissonDisc([n,n],sp);
nb = size(X,1);

F = zeros(n);
k = 30*4; % nb replication
kdisp = 10;
rpause = 5; % pause time in btw
zm = 3; % zoom on frequency
zp = 4; % zoom on points

FrKeep = [];
for it=1:k
    progressbar(it,k);
    % draw samples
    switch name
        case 'poisson-disk'
            X = poissonDisc([n,n],sp);
        case 'poisson'
            X = 1+rand(nb,2)*n;
        case 'square'
            m = linspace(-1.5,1.5,round(sqrt(nb)*3/2));
            [y,x] = meshgrid(m,m); x = (x+1i*y)*exp(2i*pi*rand); 
            x = x(abs(real(x))<1 & abs(imag(x))<1);
            X = [real(x(:))/2+.5, imag(x(:))/2+.5];
            X = X*(n-1)+1;
        case 'hexagonal'
            m = linspace(-2,2,round(sqrt(nb)*2));
            [y,x] = meshgrid(m,m); 
            y = y + x*cos(pi/3);            
            x = (x+1i*y)*exp(2i*pi*rand); 
            x = x(abs(real(x))<1 & abs(imag(x))<1);
            X = [real(x(:))/2+.5, imag(x(:))/2+.5];
            X = X*(n-1)+1;
    end
    % power spectrum
    Xi = floor(X);
    A = zeros(n);
    A( Xi(:,1) + (Xi(:,2)-1)*n ) = 1;
    f = fft2(A);
    f(1)=0;
    % accumulate
    F = F+abs(f).^2;
	Fr = ComputeRadialPofile(abs(f).^2);
    FrKeep(:,end+1) = Fr;
    %%% DISPLAY %%%
    if it<=kdisp
        s = (it-1)/(kdisp-1);        
        % plot samples
        clf;
        plot(2*X(:,1)/n-1, 2*X(:,2)/n-1, '.', 'Color', [s 0 1-s], 'MarkerSize', 15);
        axis equal; axis([-1 1 -1 1]/zp);
        set(gca, 'XTick', [], 'YTick', []); box on;
        drawnow;
        for jt=1:rpause
            saveas(gcf, [rep 'samples-' znum2str((it-1)*rpause+jt, 2) '.png']);
        end
        % colormap
        m = linspace(0,1,r+1)';
        CM = m*[s 0 1-s] + (1-m)*[1 1 1];
        %
        clf;
        imagesc([-1 1],[-1 1], fftshift(abs(f)));
        axis equal; axis([-1/zm 1/zm -1/zm 1/zm]); axis off;
        colormap(CM);
        drawnow;
        for jt=1:rpause
            saveas(gcf, [rep 'spectrum-' znum2str((it-1)*rpause+jt, 2) '.png']);
        end
        % profile
        clf;
        plot(Fr(1:round(length(Fr)/zm)), 'LineWidth', 2, 'Color', [s 0 1-s]);        
        axis([1 length(Fr)/zm 0 max(Fr)]);
        axis tight;
        set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/3 1]);
        for jt=1:rpause
            saveas(gcf, [rep 'radial-' znum2str((it-1)*rpause+jt, 2) '.png']);
        end
    end
end
F = F/k;
%
% AutoCrop(rep, 'samples')
% AutoCrop(rep, 'spectrum')
% AutoCrop(rep, 'radial')
% AutoCrop(rep, 'expectation-spectrum')
% AutoCrop(rep, 'expectation-radial')

clf;
imagesc([-1 1],[-1 1], fftshift(sqrt(F)));
axis equal; axis([-1/zm 1/zm -1/zm 1/zm]); axis off;
colormap(1-gray(256));
saveas(gcf, [rep 'expectation-spectrum.png']);

Fr = ComputeRadialPofile(F);
clf;
shadedErrorBar(1:length(Fr),Fr,4*std(FrKeep,[], 2),'k');
% plot(Fr(1:round(length(Fr)/zm)), 'LineWidth', 2, 'Color', 'k');
axis tight;
axis([1 length(Fr)/zm 0 max(Fr+4*std(FrKeep,[], 2))]);
set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/3 1]);
saveas(gcf, [rep 'expectation-radial.png']);