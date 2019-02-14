%%
% Plot of Delaunay triangulation


addpath('../toolbox/');
rep = MkResRep('del-vor');


%%
% First we create a 2D closes polygon.

clf; hold on;
f0 = [];
while true
    axis equal; axis([0 1 0 1]);  % axis off;
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    f0(end+1) = a + 1i*b;
    plot(f0, 'r', 'LineWidth', 2);
end
%%
clf; hold on;
k = size(f0,2);
plot(f0([1:end,1]), 'r', 'LineWidth', 2);
f1 = [];
for it=1:k
    axis([0 1 0 1]); axis equal; axis off;
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    f1(end+1) = a + 1i * b;
end
%%
clf; hold on;
k = size(f0,2);
plot(f0([1:end,1]), 'r', 'LineWidth', 2);
f2 = [];
for it=1:k
    axis([0 1 0 1]); axis equal; axis off;
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    f2(end+1) = a + 1i * b;
end
f0 = f0(:); f1 = f1(:); f2 = f2(:);



% add corners
% f0(end+1:end+4) = [1 1 0 0]'+1i*[1 0 0 1]';
% f1(end+1:end+4) = [1 1 0 0]'+1i*[1 0 0 1]';


t = linspace(0,1,100);
c = exp(2i*pi*t);

q = 75; 
% animate
for it=1:q
    s = (it-1)/(q-1);
    if s<1/3
        s1 = 3*s;
        Z = (1-s1)*f0 + s1*f1;
    elseif s<2/3
        s1 = 3*(s-1/3);
        Z = (1-s1)*f1 + s1*f2;
    else
        s1 = 3*(s-2/3);
        Z = (1-s1)*f2 + s1*f0; 
    end
    
    %
    T = delaunay(real(Z), imag(Z));
    [VX,VY] = voronoi(real(Z), imag(Z));
    V = VX + 1i*VY;
    % circum centers
    DT = delaunayTriangulation([real(Z), imag(Z)]);
    W = circumcenter(DT);
    Wc = W(:,1) + 1i*W(:,2);
    % circum radius
    r = abs(Wc - Z(T(:,1)));
    
    % Display Delaunay
    clf; hold on;
    % circum circles
    if size(Z,2)<=10
        plot(Wc, 'b.',  'MarkerSize', 30);
        for i=1:length(W)
            % plot( Wc(i)+r(i)*c, 'Color', [1 1 1]*.5, 'LineWidth', 1 );
        end
    end
    triplot( T, real(Z), imag(Z), 'r.-', 'LineWidth', 2, 'MarkerSize', 25 );
    axis equal; axis([0 1 0 1]);
    axis on; box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'del-'  znum2str(it,2) '.png']);

    
    % Display Delaunay
    clf; hold on;
    plot(V, 'b-', 'LineWidth', 2);
    plot(Z, 'r.',  'MarkerSize', 30);
    axis equal; axis([0 1 0 1]); 
    axis on; box on;
    set(gca, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'vor-'  znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'vor');
