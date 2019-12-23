%%
% Game of life

n = 64;

X = zeros(n);

if 1
X(end/2-1:end/2+1,end/2-1:end/2+1) = randn(3)>0;
X = randn(n)>0;
else
    X(end/2,end/2)=1;
    X = imread('CA02.jpg');
    X = (X(:,:,1)>128);
    n = size(X,1);
end


clf;
X = zeros(n);


% click selection
if 0
clf; hold on;
while true
    axis equal; axis([1 n 1 n]);
    imagesc(X'); axis equal; axis off;
    colormap gray;
    box on;
    [a,b,button] = ginput(1);
    a = round(a); b = round(b);
    if button==3
        break;
    end
    X(a,b) = 1;
end
end

X = rand(n)>.8;

s = [2:n 1];
t = [n 1:n-1];

q = 90;
Xdraw = zeros(n);
for it=1:q
 	% display
    Xdraw = X + .3*Xdraw;
    clf; 
    imagesc(1-Xdraw'); axis tight;  axis equal;  axis off;
    drawnow; 
    colormap gray(256);
    % save
    imwrite(rescale(Upsc(rescale(1-Xdraw'),2)), [rep 'anim-' znum2str(it,3) '.png']);
    % count #neighbor
    M = X(s,s) + X(s,:) + X(s,t) + ...
        X(:,s) + X(:,t) + ...
        X(t,s) + X(t,:) + X(t,t);
    X = ( (X==1) & ( M==2 | M==3 ) ) | ( (X==0) &  M==3  );
end