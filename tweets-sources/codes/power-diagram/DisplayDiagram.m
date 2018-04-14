function DisplayDiagram(Y,w,wdisp)

k = size(Y,2);
n = 400; % size of the image
N = n*n;


% coloring
C = distinguishable_colors(k+1);
C(4,:) = [];

t = linspace(0,1,n);
[v,u] = meshgrid(t,t);
X = [u(:)';v(:)'];

% weighted voronoi cells
D = distmat(X,Y).^2 + repmat(w(:)', [N 1]);
[D1,I] = min(D,[], 2);
% display

J = reshape(I, [n n]);
D1 = reshape(D1, [n n]);
clf; hold on;

% render in colors
R = zeros(n,n,3);
for m=1:k
    for ic=1:3
        R(:,:,ic) = R(:,:,ic) + (J==m)*C(m,ic);
    end
end
s = .5;
imagesc(t,t,s + (1-s)*R); colormap jet(256);

warning off;
for l=1:k
    contour(t,t,J==l,[.5 .5], 'k', 'LineWidth', 2);
end
warning on;

%for m=1:k
%    plot(Y(2,m), Y(1,m), '.', 'MarkerSize', 25, 'Color', .8*C(m,:));
%end


s = wdisp*400; % size
scatter( Y(2,:), Y(1,:), s, .8*C, 'filled' );


axis equal;
plot([0 1 1 0 0], [0 0 1 1 0], 'k', 'LineWidth', 2);
axis off;  drawnow;

end