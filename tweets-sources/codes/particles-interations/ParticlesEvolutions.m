%%
% Evolves particles according to interaction laws.

addpath('../toolbox/');
rep = MkResRep();

n = 40;
% initial positions


X = rand(n,1) + 1i*rand(n,1);
V = zeros(n,1);
V = randn(n,1)+1i*randn(n,1);

X = [-2+rand(n/2,1) + 1i*rand(n/2,1); ...
     +2+rand(n/2,1) + 1i*rand(n/2,1)]; 
X = rand(n,1) + 1i*rand(n,1);
V = zeros(n,1);
V = 0*(randn(n,1)+1i*randn(n,1));

kappa = 2; 
kappa = 3; % Coulomb/gravitation


dt = .01;
sf = 15; % skeep frame
dt = .01/sf;
A = 1; % BB

G = 1; % gravitation

rho = .01; % to remove singularity
% confinement force
M = 15*2;

PerBound = @(z)mod(real(z),1) + 1i*mod(imag(z),1);
mymod  = @(z)mod(z+1/2,1)-1/2;
mymodC = @(z)mymod(real(z)) + 1i*mymod(imag(z));


X = X-mean(X);

Col = distinguishable_colors(n);

disp_mode = 3; % 1=fast, 2=medium, 3=slow

m = 10*sf; % keep previous steps
Xold = repmat(X, [1 m]); 
it = 0;
nb = 0;
while true
    it = it+1;
    % display
    if mod(it,sf)==1
        nb = nb+1;
        if nb>300
            break;
        end
        clf; hold on;
        switch disp_mode
            case 1
                plot(X, 'k.', 'MarkerSize',25);
                plot(permute(Xold,[2,1]), 'k', 'LineWidth', 1);
            case 2
                for s=1:n
                    plot(X(s), '.', 'color', Col(s,:), 'MarkerSize',25);
                    plot(Xold(s,:), 'color', Col(s,:), 'LineWidth', 1);
                end
            case 3
                for s=1:n
                    plot(X(s), '.', 'color', Col(s,:), 'MarkerSize',25);
                    for r=1:m-1
                        u = (r-1)/(m-1);
                        plot(Xold(s,r:r+1), 'color', (1-u)*Col(s,:)+u, 'LineWidth', 1);
                    end
                end
                
                
        end
        axis equal;
        axis([-1 1 -1 1]*A); axis on; box on;
        set(gca, 'Xtick', [], 'Ytick', []);
        drawnow;
        saveas(gcf, [rep 'simul-' znum2str(nb,3),'.png']);
    end
    % evolves
    XX = repmat(X,[1,n]); 
    Delta = ( XX-permute(XX,[2 1]) );
    F = -G * sum(   Delta ./ ( (rho + abs(Delta)).^kappa ), 2 );
    F = F - M*X; % confinement force
    [X,V] = deal(X + dt * V, V + dt*F);    
    % trajectory
    Xold = [X Xold];
    Xold = Xold(:,1:m);
end

% AutoCrop(rep, 'simul-');