%%
% draw triangles as spirals

if not(exist('test'))
    test = 0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));


% input points
clf; hold on;
Z = [];
while true
    axis equal; axis([0 1 0 1]); 
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(:,end+1) = [a;b];
end
% add corner
Z(:,end+1:end+4) = [[1 1 0 0];[1 0 0 1]];
Zc = Z(1,:) + 1i*Z(2,:);

% compute delaunay
T = delaunay(Z(1,:),Z(2,:));
V = reshape(Z(1,T) + 1i*Z(2,T), size(T));


col = 'r';
lw = 1;

q = 50;
delta_list = linspace(.5,0.01,q);
delta_list = rescale(linspace(1,0,q).^3,.015,1);

for i=1:size(V,1)
    V(i,:) = V(i,randperm(3));
end

if 0 % reverse orientation
for i=1:size(V,1)
    if mod(i,2)==1
        V(i,:) = V(i,end:-1:1);
    end
end
end

for it=1:q
    t = (it-1)/(q-1);
    delta = delta_list(it);
    clf;
    DrawTriangle(V, delta, [t 0 1-t], lw);
    drawnow;
    saveas(gcf, [rep znum2str( it,2 ) '.png']);
end

% AutoCrop(rep, '');