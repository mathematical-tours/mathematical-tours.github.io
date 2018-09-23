% Test for quadree image approximation using constant 
% approximation on squares.

addpath('../toolbox/');
addpath('./toolbox-cart/');
rep = MkResRep();

name = 'hibiscus';
n = 512;
f = load_image(name, n);
f = rescale(sum(f,3));

% lagrangian penality for each node in the tree
t = 1;  % 1000
% number of scales 
J = log2(n)-1;

q = 1000;
tlist = [linspace(100,0,q-2)];
q = length(tlist);

hash_prev = 0; Z = randn(n);
idisp = 0;
for it=1:q
    t = tlist(it);    
    % error function for constant approximation by the mean
    err = @(P,d)sum( (P-repmat(mean(P,3), [1 1 size(P,3)])).^2, 3 ) + t^2;
    % generate tree 
    T = build_quadtree(f, J, err);
    W = perform_cart(T);    
    Approx = @(P,d)repmat(mean(mean(P,1),2), [size(P,1) size(P,2) 1]);
    f1 = perform_quadtree_application(f,W,Approx);
    % hash function 
    hash_new = sum(f1(:).*Z(:));
    if hash_new~=hash_prev
        idisp = idisp+1;
        % Display as a tree.
        clf; plot_tree(W);
        saveas(gcf, [rep 'tree-' znum2str(idisp,3) '.png']);
        % Display as a segmentation.
        clf;
        plot_quadtree(W,f1);
        colormap gray(256);
        caxis([0 1]);
        drawnow;
        saveas(gcf, [rep 'approx-' znum2str(idisp,3) '.png']);
    end
    hash_prev = hash_new;
end
% AutoCrop(rep, 'approx-');
% AutoCrop(rep, 'tree-', 50);