%%
% Display displacement interpolation of dense and discrete measures.

addpath('../toolbox/');
addpath('./mexEMD/');
addpath('./toolbox-lsap/');
addpath('./img/');

names = {'disk', 'two_disks'};
names = {'disk', 'joined-balls'};
names = {'cat', 'annulus'};

str = [names{1} '-' names{2}];
rep = MkResRep(str);




% helpers
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
ms = 15;
myplot = @(x,y,col)plot(x,y, '.', 'MarkerSize', ms, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'LineWidth', 2);
myplotS = @(x,y,ms,col)plot(x,y, '.', 'MarkerSize', ms, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'LineWidth', 1);


% load the shapes
n0 = 140;
n0 = 100;
xy = {};
img = {};
for i=1:2
    f = load_image(names{i},n0);
    f = sum(f,3);
    f = rescale(f)>.5;
    if f(1)==1
        f=1-f;
    end
    img{i} = f;
    I = find(f(:));
    [x,y] = ind2sub([n0 n0],I);
    xy{i} = [x(:) y(:)];
    xy{i} = (xy{i}-1)/(n0-1);
end
N = [size(xy{1},1), size(xy{2},1)];

% use fast LSAP code
use_lsap = 0;

% target points
if use_lsap
    Ntgt = min(5000,min(N(:)));
    for i=1:2
        % subsample
        I = randperm(size(xy{i},1));
        xy{i} = xy{i}(I(1:Ntgt),:);
    end
    N = [size(xy{1},1), size(xy{2},1)];
end


%%
% Display the point clouds.

clf;
subplot(2,1,1);
myplot(xy{1}(:,1), xy{1}(:,2), 'r');
subplot(2,1,2);
myplot(xy{2}(:,1), xy{2}(:,2), 'b');



q = 80; % #frames
tlist = linspace(0,1,q);
if 1
    % just linear interpolation
    for k=1:length(tlist)
        t=tlist(k);
        col = [1-t;0;t];
        B = (1-t)*img{1} + t*img{2}; 
        Br = [];
        for r=1:3
            Br(:,:,r) = B * col(r) + (1-B);
        end
        clf; imageplot(Br);
        imwrite(rescale(Br), [rep 'interp-dens-' znum2str(k,2) '.png']);
        drawnow;
    end
    return;
end



%%
% Solve the linprog of OT between the dense clouds.

meth = 'flow-matching';
meth = 'ot-lsap';
meth = 'ot';

a = ones(N(1),1) / N(1);
b = ones(N(2),1) / N(2);
C = distmat(xy{1}',xy{2}').^2;
switch meth
    case 'ot-lsap'
        C1 = int32( round(C*1e6) );
        [J,varrho,u,v] = hungarianLSAP(C1);
        I = (1:N(1))';
        gammaij = ones(N(1),1)/N(1);
    case 'ot'
        [cost,gamma] = mexEMD(a,b,C);
        [I,J,gammaij] = find(gamma);
    case 'flow-matching'
        [I,J] = meshgrid(1:N(1),1:N(2));
        I = I(:); J = J(:);
        Ntgt = round(.2*prod(N)); % number of keps points
        s = randperm(prod(N)); s = s(1:Ntgt); % sub-sample
        I = I(s); J = J(s);
        gammaij = ones(prod(Ntgt),1)/Ntgt;
end

%%
% Render using density estimator.


for k=1:length(tlist)
    t=tlist(k);
    col = [1-t;0;t];
    Xt = (1-t)*xy{1}(I,:) + t*xy{2}(J,:);
    %%% render as point clouds
    if 0
        Nt = length(I);
        s = ones(Nt,1)*20; % size
        col = cat(2,(xy{1}(I,:)), ones(Nt,1));
        % col = cat(2,cos(10*pi*xy{1}(I,:)), ones(Nt,1));

        u = xy{1}(I,1); v = xy{1}(I,2);
        m = 6;

        t = (mod(u*m,1)<.5 ) == (mod(v*m,1)<.5);
        t = rescale( cos(10*pi*u)  .* sin(10*pi*v) );

        t = xy{2}(J,1);

        col = (1-t)*[1 0 0] + t*[0 0 1];

        clf; scatter( Xt(:,2), Xt(:,1), s, col, 'filled' );
        axis equal; axis([0 1 0 1]); box on; set(gca, 'XTick', [], 'YTick', []);
        axis ij;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(k,2) '.png']);
    else
        %%% render as image
        K = 10; % up-sampling
        n1 = n0*K;
        Xt1 = round( K*Xt*(n0-1) + 1);
        A = zeros(n1,n1);
        %        A(Xt1(:,1) + n1*(Xt1(:,2)-1)) = gammaij;
        for m=1:size(Xt1,1)
            A(Xt1(m,1) + n1*(Xt1(m,2)-1)) = A(Xt1(m,1) + n1*(Xt1(m,2)-1)) + ...
                gammaij(m);
        end
        nF = 128*2; % width of the convolution kernel
        s = K*2;
        s = K*.75;
        g = exp(-(-nF:nF).^2 / (2*s^2) );
        B = conv2(conv2(A, g, 'same')', g, 'same')';
        %
        Br = [];
        for r=1:3
            Br(:,:,r) = B * col(r) + (1-B);
        end
        clf; imageplot(Br);
        imwrite(rescale(Br), [rep 'interp-dens-' znum2str(k,2) '.png']);
        drawnow;
    end
end

% AutoCrop(rep, 'anim');

return;

%%
% Render on a small cloud of discrete samples.

P = 300;
for i=1:2
    I = randperm(N(i)); I = I(1:P);
    xyS{i} = xy{i}(I,:);
end

CS = distmat(xyS{1}',xyS{2}').^2;
[cost,gammaS] = mexEMD(ones(P,1)/P,ones(P,1)/P,CS);

[I,J,gammaij] = find(gammaS);
tlist = linspace(0,1,q);
ms = 30;
for k=1:length(tlist)
    t=tlist(k);
    Xt = (1-t)*xyS{1}(I,:) + t*xyS{2}(J,:);
    % plot as points
    clf;
    hold on;
    for i=1:length(gammaij)
        myplotS(Xt(i,2), Xt(i,1), ms, [1-t 0 t]);
    end
    axis([0 1 0 1]);
    % dummy points
    plot([0 1], [0 1], '.', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
    axis equal; axis square;
    axis off; axis ij;
    drawnow;
	saveas(gcf, [rep 'interp-points-' znum2str(k,2) '.png'], 'png');
end

% AutoCrop(rep, 'interp-points-');
