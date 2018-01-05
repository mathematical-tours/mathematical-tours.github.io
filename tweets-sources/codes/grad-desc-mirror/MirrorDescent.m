%%
% Mirror descent

rep = '../results/mirror-descent/';
[~,~] = mkdir(rep);

bregmode = 'kl';
bregmode = 'sqrt';
bregmode = 'log';
bregmode = 'eucl';


c = [1 4]';

fname = 'linear';
fname = 'quadratic';

switch fname
    case 'quadratic'
        f = @(u,v)(u.^2+v.^2)/2;
        Gradf = @(x)x;
    case 'linear'
        f = @(u,v)u*c(1)+v*c(2);
        Gradf = @(x)c;
end

N = 256;
tx = linspace(0,1,N);
ty = linspace(0,.6,N);
[Y,X] = meshgrid(ty,tx);
F = f(X,Y);

clf; hold on;
imagesc(tx,ty,-F');
contour(tx,ty,F', 15, 'k');

col = {'r' 'g' 'k' 'y'};
breglist = {'kl' 'sqrt' 'log' 'eucl'};

tau = .05;
niter = round(100/tau);

for it = 1:length(bregmode)

    bregmode = breglist{it};

    [R,Ri,phi] = load_bregman(bregmode);


    x = [.95 .55]';
    Z = [];
    for i=1:niter
        Z(:,end+1) = x;
        x = Ri( R(x) - tau*Gradf(x) );
        x = max(x,0);
    end


    plot(Z(1,:),Z(2,:), '-', 'LineWidth', 2, 'MarkerSize', 23, 'color', col{it});
    axis off; axis equal;

end


saveas(gcf, [rep 'mirror-' fname '.png'], 'png');



return;

% triangle simplex

W = 1-U-V;

normalize = @(p)p/sum(p);
p = normalize([1 3 1]);

KL = @(X,Y) sum( X.*log(X./Y) + Y - X, 2);

a = KL( [U(:) V(:) W(:)], repmat(p(:)', [n*n 1]) );
a = reshape(a, [n n]);
a(W<=0)=NaN;

b = KL( repmat(p(:)', [n*n 1]), [U(:) V(:) W(:)] );
b = reshape(b, [n n]);
b(W<=0)=NaN;
