%%
% integrate a vector field.


s = 100;
s = 20;
s = 40;

addpath('../toolbox/');
rep = MkResRep(num2str(s));


normVF = @(vx)vx ./ repmat( sqrt(sum(vx.^2,2)),[1 2] );

n = 400; 

% gaussian blur
t = [0:n/2,-n/2+1:-1]';
normalize = @(x)x/sum(x(:));
G = @(s)normalize(exp(-t.^2/(2*s^2))); G = @(s)G(s)*G(s)';
fconv = @(x,y)real( ifft2( fft2(x).*fft2(y) ) );
GFilt = @(f,s)fconv(f, G(s));

% generate a random vector field.
v = cat(3, GFilt(randn(n),s), GFilt(randn(n),s));

a = sqrt(sum(v,3));
v = v/std(a(:));

% add dummy row
v(:,end+1,:) = v(:,1,:);
v(end+1,:,:) = v(1,:,:);

t = linspace(0,1,n+1);
[Yi,Xi] = meshgrid(t,t);


% seed points
p = 200; 
x0 = rand(p,2);

p = 15;
t = linspace(0,1,p+2); t = t(2:end-1);
[Yg,Xg] = meshgrid(t,t);
x0 = [Xg(:),Yg(:)];

% for the display of the vectors
% p = 15*2;
t = linspace(0,1,p+2); t = t(2:end-1);
[Yg,Xg] = meshgrid(t,t);
v1 = InterpVF(Xi,Yi,v,[Xg(:),Yg(:)]);
v1 = normVF(v1);
v1 = reshape(v1, [p p 2]);

% double (2 directions)
x0 = [x0;x0];


T = .2;
dt = .01/10;
niter = round(T/dt);

q = 50;
ndisp = round(linspace(1,niter,q));
kdisp = 1;

x = x0;
X = [x(:,1)+1i*x(:,2)];
for it=1:niter
    
    % interpolate positions
    vx = InterpVF(Xi,Yi,v,x);
    % normalize
    vx = normVF(vx);
    
    % update points
    x(1:end/2,:) = x(1:end/2,:) - dt*vx(1:end/2,:);
    x(end/2+1:end,:) = x(end/2+1:end,:) + dt*vx(end/2+1:end,:);
    x = mod(x,1);
    X = [X, x(:,1)+1i*x(:,2)];
    % draw;
    if it==ndisp(kdisp)
        clf; hold on;
        %
        PlotFlow(X);
        %
        quiver(Xg,Yg,v1(:,:,1), v1(:,:,2), .4, 'k', 'LineWidth', 2);
        %
        axis equal; axis([0 1 0 1]); axis off;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(kdisp,2) '.png']);
        kdisp = kdisp+1;
    end
end

% AutoCrop(rep, 'anim');

