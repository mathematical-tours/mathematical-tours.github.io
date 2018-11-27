%% 
% implement the Cole-Hopf transformation h=exp(f/epsilon) 
%      |nabla f|^2 - epsilon*(Delta f) = 1	(*)
% -->  Delta h = -h/epsilon

n = 150;

name = 'elephant';
name = 'triskel';
name = 'dots';
name = 'cat';

addpath('../toolbox/');
rep = MkResRep(name);

if strcmp(name, 'dots')
    rho = ones(n);    
    rand('state', 1342);
    S = round( rand(6,2)*n );
    eps_range = [.14 10000];
    c = [];
else
	f = load_image(name, n);
    f = sum(f,3); f = double((f/max(f(:)))>.5);
    if f(1)==0
        f = 1-f;
    end
    % click selection
    if not(exist('S'))
    S = [];
    clf; imagesc(f); hold on;
    while true
        axis equal; axis([1 n 1 n]);
        box on; set(gca, 'Xtick', [], 'Ytick', []);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if button==3
            break;
        end
        S(end+1,:) = round([b a]);
    end
    end
    % should be small inside
    rho = (1-f)*1e5 + f;
    [C,h] = contour(f,[.5,.5]);
    m = C(2,1);  c = C(:,2:m+1); % c = (c-1)/(n0-1);
    %    
    eps_range = [.05 10000];
    eps_range = [.05 50];
end


% grad operator with neumann BC
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);
dx = dx(2:end,:);

% with periodic BC
dx = spdiags([-ones(n,1), ones(n,1)], -1:0, n, n);
dx(1,end) = -1;

%
Grad = [kron(dx,speye(n)); kron(speye(n),dx)];
Delta = Grad'*Grad / n^2;

Delta = Grad' * spdiags( [rho(:); rho(:)], 0, 2*n*n,2*n*n ) * Grad / n^2;


% set of points and diracs RHS
I = S(:,1)+(S(:,2)-1)*n;
b = zeros(n*n,1); b(I) = 1; % dirac masses


q = 50;
eps_list = eps_range(1) + (eps_range(2)-eps_range(1))*linspace(0,1,q).^4;

for it=1:q
    s = (it-1)/(q-1);
    epsilon = eps_list(it);
    %
    U = Delta+speye(n*n)/epsilon;
    U(I,:) = 0; U(I+(I-1)*n*n) = 1; % eye on diag    
    % solve
    h = U\b;
    h = reshape(h,[n n]);
    D = -epsilon*log(h);
    if not(isempty(c))
        D(f==1) = NaN;
        D = D/max(D(f==0));
        %
        r = 20; % #levellines
    else
        D = rescale(D);
        r = 15; % #levellines
    end
    m = linspace(0,1,r-1)';
    CM = m*[s 0 1-s] + (1-m)*[1 1 1];
    %
    clf; hold on;
    imagesc(D);
    contour(D,linspace(0,1,r), 'Color', [s 0 1-s]*.8, 'LineWidth', 2);
    colormap(CM);
    caxis([0 1]);
    axis image; axis off;
    plot(S(:,2),S(:,1), '.r', 'MarkerSize', 25);
    axis image; axis off;
    if not(isempty(c))
        plot(c(1,:), c(2,:), 'color', 'k', 'LineWidth', 2);
    end
    axis ij;
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    
end

% AutoCrop(rep,'anim');

