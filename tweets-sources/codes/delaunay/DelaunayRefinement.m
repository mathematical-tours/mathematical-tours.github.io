%%
% Delaunay refinement aka Ruppert's algorithm

rep = '../results/delaunay-refinement/';
[~,~] = mkdir(rep);
addpath('../toolbox/');


% input points
clf; hold on;
Z = [];
Z = [[1 1 0 0];[1 0 0 1]];
while true
    axis([0 1 0 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(:,end+1) = [a;b];
end
% add corner
Zc = Z(1,:) + 1i*Z(2,:);

niter = 50;

t = linspace(0,1,100);
c = exp(2i*pi*t);

% constraints
Constr = [[1 2]; [2 3]; [3 4]; [4 1]];

e = 0;
Proj = @(x)min(max(x,+e),1-e);

for i=1:niter
	% circum centers
    warning off;
    DT = delaunayTriangulation([real(Zc);imag(Zc)]', Constr);
    warning on;
    T = DT.ConnectivityList;
    V = DT.Points(:,1) + 1i*DT.Points(:,2);
    W = circumcenter(DT);
    Wc = W(:,1)' + 1i*W(:,2)';
    % circum radius
    r = abs(Wc - Zc(T(:,1)));
    % only those inside
    if 0
    I = find(real(Wc)<1+e & real(Wc)>-e & imag(Wc)<1+e & imag(Wc)>-e);
    Wc = Wc(I); r = r(I);
    end
    % farthest
    [~,j] = max(r);
    w = Wc(j);
	w = Proj(real(w))+1i*Proj(imag(w)); 
    
    
    % Display 
	clf; hold on;
    triplot( T, real(Zc), imag(Zc), 'r.-', 'LineWidth', 2, 'MarkerSize', 25 );
    axis equal; axis([0 1 0 1]);
    axis off;
    % circum circles
    plot(Wc(j), 'k.', 'MarkerSize', 25);
    plot(w, 'k.', 'MarkerSize', 25);
    plot( Wc(j)+r(j)*c, 'k', 'LineWidth', 2 );
    saveas(gcf, [rep 'ruppert-' num2str(i) ' .eps'], 'epsc');
    
    % add farthest
    Zc(end+1) = w;

end

