

name = 'france';
ext = 'jpg';

name = 'triskel';
name = 'elephant';
name = 'horse';
ext = 'png';

addpath('../toolbox/');
rep = MkResRep(name);


% Generate random planar graph

n = 1000; 

f = sum(imread([name '.' ext]),3);
f = f>mean(f(:));
if f(1)==1
    f = 1-f;
end

[X,Y] = meshgrid(linspace(0,1,size(f,2)), linspace(0,1,size(f,1)));
I = find(f==1); I = I(randperm(length(I)));
P = [X(I(1:n)) Y(I(1:n))];


%
DT = delaunayTriangulation(P);
T = DT.ConnectivityList;
I = T(:);
J = [T(:,2); T(:,3); T(:,1)];
d = ones(length(I),1); % combinatorial
d = sqrt( sum( (P(I,:) - P(J,:)).^2,2 ) ); % distances
A = sparse(I,J,d,n,n);
A = max(A,A');

% A(A>.06) = 0; % elephant
A(A>.06) = 0; % triskel

clf; hold on;
gplot(A,P, 'k');
plot(P(:,1), P(:,2), 'b.', 'MarkerSize', 25);
axis equal; axis tight; axis ij;

z = .5+.5i;
switch name
    case 'elephant'
        z = .95 + .81i;
    case 'triskel'
        z = .74  + .74i;
    case 'horse'
        z = .9 + .8i;
end

[v,Istart] = min(abs(P(:,1)+1i*P(:,2)-z)); Istart = Istart(1);
% target point
[d,Iend] = max(abs(P(:,1)+1i*P(:,2)-z)); Iend = Iend(1);
h = abs(P(:,1)+1i*P(:,2)-P(Iend,1)-1i*P(Iend,2));


clf; hold on;
gplot(A,P, 'k');
plot(P(Istart,1), P(Istart,2), 'b.', 'MarkerSize', 25);
plot(P(Iend,1), P(Iend,2), 'r.', 'MarkerSize', 25);
axis equal; axis tight; axis ij;


opt.svg_rate = 1; % ceil(n/200);
opt.target = Iend;
opt.heuristic = 2*h; % set to 0 for no heuristic
[D0,Dsvg,Ssvg] = dijkstra(A, Istart, opt);

q = 50; 
ndisp = round(linspace(1,size(Ssvg,2),q));
for it=1:q
    i = ndisp(it);
    S = Ssvg(:,i);
    D = Dsvg(:,i);
    %
    I = find(D~=Inf);
    J = find(D==Inf);
    clf; hold on;
    gplot(A,P, 'k');
    %
    s = ones(length(J),1)*60; % size
    c = ones(length(J),3)*.5;
    scatter( P(J,1), P(J,2),s,c,'filled');
    %
    s = ones(length(I),1)*80; % size
    m = D(I)/D0(Iend);
    c = m*[1 0 0] + (1-m)*[0 0 1]; % colors
    scatter( P(I,1), P(I,2),s, c, 'filled');
    %
    axis equal; axis([0 1 0 1]); box on;  axis ij;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    %
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end




% AutoCrop(rep, 'anim');