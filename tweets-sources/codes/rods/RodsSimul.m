%%
% Simulation of a system of rods.

name = 'horizontal';
name = 'rand';


addpath('../toolbox/');
rep = MkResRep(name);

n = 16; % #points


% min E(x) = 1/2 * | abs(D*x)-d |^2
% E(x+e) =  1/2 * | abs(D*x)-d + (D*e) .* (D*x)./|D*x| |^2
%   = E(x) + <abs(D*x)-d, (D*e) .* (D*x)./|D*x|>
% nabla E(x) = D' * [ (abs(D*x)-d) .* (D*x)./|D*x| ]
D = eye(n)-diag(ones(n-1,1),1);
D = D(1:n-1,:);


u = linspace(0,1,n)';


srand = @(n)2*rand(n,1)-1;

switch name
    case 'horizontal'
        I = [1];
        x = u;
        BB = [-1.1 1.1 -1.1 .1];
        
    case 'diagonal'
        I = [1];
        x = .3*u + .8*1i*u;
        BB = [-.4 1.1 -1.1 .1];
        BB = [-.4 1.1 -sum(abs(D*x))-.05 1];
        
    case 'rand1'
        x = srand(n)*.4 + 1i*srand(n)*.02;
        x(1)=-.5; x(end)=+.5;
        BB = [-1 1 -sum(abs(D*x))-.05 .1];
        
        
    case 'rand'
        x = srand(n)*.35 + 1i*srand(n)*.05;
        I = [1 round(.7*n)];
        x(I)=[-.5,+.5];
        BB = [-.8 1.3 -1.75 .1];
end


% fixed position
xI = x(I);

% gradivity forces
f = (0-1i)*10;

Tmax = 1.75;
niter = 300;
dt = Tmax/niter;

q = 80; % #frames


% for inner optim
tau = 5;
niter_proj = 1000;

ndisp = round(linspace(1,niter,q));
k = 1;

% evolution
v = zeros(n,1);
d = abs(D*x); % fixed length
for it=1:niter
    % acceletate
    E = [];
    % tentative point
    z = x-dt*v;
    for j=1:niter_proj
        Dx = D*z;
        E(j) = norm(abs(Dx)-d);
        nablaE = 1/n * D' * ( (abs(Dx)-d) .* (Dx)./abs(Dx) );
        z = z - tau*nablaE;
        z(I)=xI;
    end
    % corrected speed
    v = (x-z)/dt; 
    % update
    x = x-dt*v; x(I)=xI;
    v = v-dt*f;
    
    clf; 
    DrawRod(x);
    % axis off; 
    set(gca, 'XTick', [], 'YTick', []); box on;
    axis equal;
    axis(BB);  
    drawnow;
    % save image
    if it==ndisp(k)
        saveas(gcf, [rep 'img-' znum2str(k,2) '.png']);
        k = k+1;
    end
    
end

% AutoCrop(rep, 'img');