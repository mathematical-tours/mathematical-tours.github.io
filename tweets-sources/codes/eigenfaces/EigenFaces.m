if 0
    % load data
    M = [];
    for i=1:40
        for j=1:10
            a = imread(['data/s' num2str(i) '/' num2str(j) '.pgm']);
            M(:,:,end+1) = double(a)/255;
        end
    end
save('database-att.mat', 'M');
end


addpath('../toolbox/');
rep = MkResRep();

load database-att.mat;

imwrite(mean(M,3), [rep 'mean.png']);

p = size(M,1)*size(M,2);
n = size(M,3);
Ma = reshape(M, [p n]);
Mmean = mean(Ma,2);
Mm = Ma - repmat(Mmean, [1 n]);

resh = @(f)reshape(f, [size(M,1) size(M,2)]);;

[U,S,V] = svd(Mm,'econ');
S = diag(S);

clf;
plot(S, 'k', 'LineWidth', 2);
axis tight;
set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20);
saveas(gcf, [rep 'eigenvalues.eps'], 'epsc');

K = 15;
for i=1:K
    u = resh(U(:,i));
    r = 2;
    u = (clamp(u/std(u(:)), -r,r)+r)/(2*r);
    imwrite(u, [rep 'eigen-' znum2str(i,2) '.png']);
end

q = 50; 
ndisp = round(1+ (n-1)*linspace(0,1,q).^2);

% Progressive de-compression
for it=1:4 % guy to test
    f = Mm(:,1+(it-1)*10); 
    f = f(:);    
    for i=1:q
        r = ndisp(i);
        f1 = U(:,1:r)*( U(:,1:r)'*f ) + Mmean;
        f1 = resh(f1);
        clf; imagesc(f1);
        drawnow;
        imwrite(f1, [rep 'reconstr-' znum2str(it,1) '-' znum2str(i,2) '.png']);
    end
end

for it=1:8 % guy to test
    f = Ma(:,1+(it-1)*10); 
    imwrite(resh(f), [rep 'input-' znum2str(it,2) '.png']);
end

% on the curve
for i=1:q
    r = ndisp(i);
    s = (i-1)/(q-1);
    clf; hold on;
    plot(S, 'k', 'LineWidth', 2);
    plot(r, S(r), '.', 'Color', [s 0 1-s], 'MarkerSize', 45);
    axis tight; box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,2) '.png']);
end
% AutoCrop(rep, 'anim');
