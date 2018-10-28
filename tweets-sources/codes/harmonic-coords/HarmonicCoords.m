% Generalized Barycentric Coordinates for Warping
%%

if not(exist('test'))
    test = 0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));


%%
% Background texture 

% grid
n = 200;
x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);


% texture
m  = 30;
[XT,YT] = meshgrid((0:n-1)/n,(0:n-1)/n);
T = mod( floor(XT*m)+floor(YT*m),2 );


name = 'lisa';
T = load_image(name,n);
T = clamp(T/255,0,1);

%%
% First we create a 2D closes polygon, which will be a cage used to perform
% 2D shape deformation.

clf; hold on;
imagesc(x,x,T);
V = [];
while true
    axis([0 1 0 1]); axis equal; axis off; axis ij;
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    V(:,end+1) = [a;b];
    plot(V(1,:), V(2,:), 'r', 'LineWidth', 2);
end
k = size(V,2);
plot(V(1,[1:end,1]), V(2,[1:end,1]), 'r', 'LineWidth', 2);
V2 = [];
for it=1:k
    axis([0 1 0 1]); axis equal; axis off;
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    V2(:,end+1) = [a;b];
end

%%
% Indicator of the shape.


S = 1 - inpolygon(X,Y,V(1,:),V(2,:));

%%
% Set to black outside. 

% T(S==1) = 1;

%%
% Display it.

lw = 3; ms = 25;
clf; hold on;
plot_surf_texture(cat(3,X,Y,zeros(n)), T);
h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
view(2); axis('off'); axis('equal'); axis ij;

%%
% Useful operator.

dotp = @(a,b)sum(a.*b);
crossp = @(a,b)a(1,:).*b(2,:)-a(2,:).*b(1,:);
normalize = @(a)a./repmat(sqrt(sum(a.^2)), [2 1]);

%%
% Points of the domain.

W = [X(:)';Y(:)'];

%%
% Compute the full set of mean coordinates.
C = zeros(n,n,k);
for i=1:k
    vi = V(:,i);
    U = repmat(vi,[1 n^2])-W;
    nb = normalize( U );
    % length
    d = sqrt( sum(U.^2) );
    s = 1;
    for j=mod([i-2,i],k)+1
        vj = V(:,j);
        na = normalize( repmat(vj,[1 n^2])-W );
        % sign
        si = s*sign(crossp(na,nb));
        % angle
        dp = dotp(na,nb);
        theta = si .* acos(clamp(dp,-1,1));
        % add tangent of half angle
        C(:,:,i) = C(:,:,i) + reshape( tan(theta/2) ./ d, [n n]);
        s = -s;
    end
end
% Normalize them.
C = C ./ repmat( sum(C,3), [1 1 k] );

%%
% display the functions
    
r = 12;
for i=1:k
    c = abs(C(:,:,i)) + 1e-3;
    c(S==1) = 0;
    opt.cm = parula(r-1);
    B = display_shape_function(c',opt);
    clf; hold on;
    t = linspace(0,1,n);
    imagesc(t,t,B); axis('image'); axis('off');
    contour(t,t,c', linspace(0,1,r), 'k', 'LineWidth', 2); axis ij;
    h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
    set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
    drawnow;
    saveas(gcf, [rep 'coords-meanval-' num2str(i)], 'png');
end
% AutoCrop(rep, 'coords-meanval-');


%%
% Interpolation function.

applyinterp = @(C,x)sum(repmat(reshape(x(:), [1 1 k]),[n n 1]).*C,3);

%%
% Apply the interpolation weight to the X/Y coordinate to test for the
% linear precision of the coordinates.

clf;
imageplot(applyinterp(C,V(1,:)), 'Should be X', 1,2,1);
imageplot(applyinterp(C,V(2,:)), 'Should be Y', 1,2,2);
colormap jet(256);

%%
% Modify the position of the cage.

q = 50; % #frames
for it=1:q
    t = (it-1)/(q-1);
    V1 = V*(1-t) + V2*t;
    % Warp the grid.
    X1 = applyinterp(C,V1(1,:));
    Y1 = applyinterp(C,V1(2,:));
    % Display the warped texture.
    clf; hold on;
    plot_surf_texture(cat(3,X1,Y1,zeros(n)), T);
    h = plot(V1(1,[1:end 1]), V1(2,[1:end 1]), 'r.-');
    set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
    view(2); axis('off'); axis('equal'); axis ij;
    axis([0,1,0,1]); drawnow;
    %
    saveas(gcf, [rep 'warp-' znum2str(it,2) '.png'], 'png');
end

% AutoCrop(rep, 'warp-');




