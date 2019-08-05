% test for traveling saleseman interface.

addpath('../toolbox/');
rep = MkResRep();

xy = rand(5,2);
s = .1; t = .9;
xy = [s s;s t;t s;t t]; 
q = 50;
niter = 1000;
for it=1:q
    userConfig = struct('xy',xy, 'SHOWPROG', false, 'SHOWRESULT', false, 'SHOWWAITBAR', false, 'NUMITER', niter);
    resultStruct = tsp_ga(userConfig);
    I = resultStruct.optRoute([1:end 1]);
    % 
    clf; hold on;
    plot(xy(I,1), xy(I,2), 'b', 'LineWidth', 2);
    plot(xy(:,1), xy(:,2), 'k.', 'MarkerSize', 20);
    axis equal; axis([0 1 0 1]);  box on;
    set(gca, 'XTick', [], 'YTick', []);
    %
    [a,b,button] = ginput(1);
    xy(end+1,:) = [a;b];
end

n = 100;
t = rand(n,1);
xy = .5 + (.3+.7*t).*[cos(2*t*2*pi),sin(2*t*2*pi)]*.47;
plot(xy(:,1), xy(:,2), '.');

n = 300;
xy = rand(n,2);

n = 20*20;
t = linspace(.05,.95,20);
[Y,X] = meshgrid(t,t);
I = randperm(n)';
xy = [X(I),Y(I)];

np_pt = round(linspace(5,n,q));

niter = 20000;
for it=1:q
    xy1 = xy(1:np_pt(it),:);
    progressbar(it,q);
    userConfig = struct('xy',xy1, 'SHOWPROG', false, 'SHOWRESULT', false, 'SHOWWAITBAR', false, 'NUMITER', niter);
    resultStruct = tsp_ga(userConfig);
    I = resultStruct.optRoute([1:end 1]);
    % 
    clf; hold on;
    plot(xy1(I,1), xy1(I,2), 'b', 'LineWidth', 2);
    plot(xy1(:,1), xy1(:,2), 'k.', 'MarkerSize', 20);
    axis equal; axis([0 1 0 1]);  box on;
    set(gca, 'XTick', [], 'YTick', []); drawnow;
    %
    saveas(gcf, [rep 'anim' '-' znum2str(it,2) '.png']);
end


AutoCrop(rep, 'anim');