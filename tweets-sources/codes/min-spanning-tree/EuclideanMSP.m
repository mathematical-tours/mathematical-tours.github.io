%%
% Test for Euclidean MSP.

addpath('../toolbox/');
addpath('./msp-toolbox/');
rep = MkResRep();

% generate points.
X = randn(2,100);
clf; PlotMST(X);
axis equal; axis tight; axis off;
saveas(gcf, [rep 'example.png']);

X = rand(2,3);
clf; PlotMST(X);
it = 0; button = 1;
while button==1
    it = it+1;
    saveas(gcf, [rep 'msp-' znum2str(it,3) '.png']);
    axis([0 1 0 1]);
    [a,b,button] = ginput(1);
    X(:,end+1) = [a;b];
    clf; PlotMST(X); 
end