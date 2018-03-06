%%
% Brenier in heat

rep = '../results/brenier-heat/';
[~,~] = mkdir(rep);

addpath('../toolbox/');
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

n = 256;
lev = 3;

name = 'hibiscus';
name = 'brenier';
name = 'benamou';
psi0 = sum( load_image(name, n), 3);

% equalize
[~,I] = sort(psi0(:)); psi0(I) = linspace(0,1,length(psi0(:)));
% smooth a bit
a = [2:n n];
b = [1 1:n-1];
sm = @(f)( f + f(a,:) + f(b,:) + f(:,a) + f(:,b) )/5;
for k=1:4
    psi0 = sm(psi0);
end
[~,I] = sort(psi0(:)); psi0(I) = linspace(0,1,length(psi0(:)));

switch lev
    case 2
        % 2 values
        psi0 = rescale(psi0, 0, 2);
    case 3
        % 3 values
        psi0 = rescale(psi0, 0, 3);
end


Delta = @(f)(4*f - f(a,:) - f(b,:) - f(:,a) - f(:,b) )/4;

niter = 150000*3;
ndisp = max(1,ceil(niter/20));
tau = 1;

% first run
clf; hold on;
psi = psi0;
for i=1:niter
    if mod(i,ndisp)==1
        imagesc(psi);
        axis image; axis off; colormap gray(256);
        drawnow;
    end
    psi = psi - tau * Delta(psi).*cos( pi*psi ).^2;
end
axis tight;


% svg run
k = 0;
E = norm(psi-psi0);
nextE = 0;
kdisp = 50;
%
clf; hold on;
psi = psi0;
for i=1:niter    
    if norm(psi-psi0)>=nextE
        k = k+1;
        % save
        imwrite(rescale(psi), [rep name '-lev' num2str(lev) '-' num2str(k) '.png']);
        % display 
        imagesc(psi);
        axis image; axis off; colormap gray(256);
        drawnow;
        % display levelsets
        r = 11; % #levellines
        clf; hold on;
        u = linspace(0,1,n);
        imagesc(u,u,psi);
        contour(u,u,psi,linspace(0,lev,r), 'k');
        colormap(parula(r-1));
        caxis([0 lev]);
        axis image; axis off; axis ij;
        saveas(gcf, [rep name '-lev' num2str(lev) '-lsets-' num2str(k) '.png'], 'png');
        % histogram   
        c = (k-1)/(kdisp-1);
        t = linspace(0,lev, 51);
        h = hist(psi(:), t);
        clf;
        bar(t,h, 'FaceColor', [c 0 1-c]); axis tight; box on; 
        SetAR(1/2);
        set(gca, 'XTick', [], 'YTick', []);
        saveas(gcf, [rep name '-lev' num2str(lev) '-hist-' num2str(k) '.png'], 'png');
        % update runtime
        nextE = norm(psi-psi0) + E/kdisp;
    end
    psi = psi - tau * Delta(psi).*cos( pi*psi ).^2;
end
axis tight;

