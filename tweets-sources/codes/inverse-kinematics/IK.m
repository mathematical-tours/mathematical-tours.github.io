%%
% Simple harm

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 4;
r = .04 + .4*rand(n,1);


n = 8;
r = .04 + .15*rand(n,1);


n = 6;
r = ones(n,1)/n;
randn('state', 124);
t0 = rand(n,1)*2*pi;

t0 = .04*randn(n,1)*2*pi - .7*pi;
t0 = .001*randn(n,1)*2*pi;


% target
y = .3+.3i;

x = forward_pass(r,t);
clf; draw_linkage(x,y);

[dt,dx] = grad_linkage(r,t,y);

Last = @(x)x(end);
f = @(t)1/2*abs(Last(forward_pass(r,t))-y).^2;
dt1 = fd_grad(f,t);

my_grad = @(t)fd_grad(f,t);
my_grad = @(t)grad_linkage(r,t,y);

t = t0;
tau = .001; % not tracking exactly
tau = .1; 
% for L-M
rho = 0.6;
niter_grad = 60;
if 1

    for k=1:niter_grad
        if 0
            % gradient descent
            t = t - tau*my_grad(t);
        else
            % Levenberg Marquardt
            JJ = diag(r) * cos(t-t') * diag(r);
            t = t - inv( JJ + rho * eye(n) ) * my_grad(t);
        end
    %
    clf; draw_linkage( forward_pass(r,t) ,y);
    axis equal; axis([-1 1 -1 1]*1.1);
    axis on; box on; set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(k);
end
end

return;


% animate points in a square with bouncing
rand('state', 125); randn('state', 1666);
k = 1; y = 2*rand(k,1)-1+1i*(2*rand(k,1)-1);
eta = .04; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 200;
for it=1:q
    % IK adaptation
    noise = (randn+1i*randn)*.0;
    f = @(t)1/2*abs(Last(forward_pass(r,t))-y).^2;
    my_grad = @(t)fd_grad(f,t);
    t = t+randn(n,1)*.1e-1; % to break ties and non smoothness
    for k=1:niter_grad 
        t = t - tau*my_grad(t);
    end
    clf; draw_linkage( forward_pass(r,t) ,y);
    axis equal; axis([-1 1 -1 1]*1.1);
    axis on; box on; set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % animate
    y = y + v;
    I = find( real(y)<-1 | real(y)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(y)<-1 | imag(y)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    mysaveas(it);
end
