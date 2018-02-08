%%
% Display for product of Gaussian rule

rep = '../results/hump-algebra/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

if not(exist('test'))
    test = 0;
end
test = test+1;

n = 256*2;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
s = .15; % width of the Gaussians
G = @(z)exp( (-(X-real(z)).^2-(Y-imag(z)).^2)/(2*s^2) );

% levelset for drawing
r = 15*2;
col = {'k', 'k', 'm', 'g', 'b'};
Z = {};
for k=1:2
    clf; hold on;
    Z{k} = [];
    M{k} = zeros(n);
    it = 0;
    while true
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if button==3
            break;
        end
        Z{k}(end+1) = a+1i*b;
        M{k} = M{k} + (-1)^it * G(Z{k}(end));
        it = it+1;
    end
    % render
    clf; hold on;
    imagesc(t,t,M{k}');
    contour(t,t,M{k}',linspace(-1,1,r), 'k');
    colormap(jet(r-1));
    caxis([-1 1]);
    plot(real(Z{k}), imag(Z{k}), '.', 'color', col{k}, 'MarkerSize', 30);
    axis tight; axis equal; axis off;
    saveas(gcf, [rep 'hump-' num2str(test) '-' num2str(k) '.png'], 'png');
end

% compute center
W = [];
for i=1:length(Z{1})
    for j=1:length(Z{2})
        W(end+1) = (Z{1}(i)+Z{2}(j))/2;
    end
end

% render
MM = M{1}.*M{2};
clf; hold on;
u = max(abs(MM(:))); 
imagesc(t,t,MM');
contour(t,t,MM',linspace(-u,u,r), 'k');
colormap(jet(r-1));
caxis([-u u]);
for k=1:2
    plot(real(Z{k}), imag(Z{k}), '.', 'color', col{k}, 'MarkerSize', 20);
end
plot(real(W), imag(W), '.', 'color', 'k', 'MarkerSize', 30);
axis tight; axis equal; axis off;
saveas(gcf, [rep 'hump-' num2str(test) '.png'], 'png');