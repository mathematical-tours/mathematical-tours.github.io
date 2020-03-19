
addpath('../toolbox/');
rep = MkResRep();

name = 'duck';
name = 'armadillo';
name = 'david_head';
name = 'gargoyle';
name = 'duck';
name = 'elephant';
name = 'bunny';
[X,F] = read_off([name '.off']);
n = size(X,2);
clear options; options.name = name;

% center
X = X - repmat( mean(X,2), [1 n] );
X = X / mean(sqrt(sum(X.^2)));

options.face_vertex_color = ones(n,1);
clf; 
plot_mesh(X,F,options);
camlight;


% Laplacian with cotan weights
lapl = 'cotan';
lapl = 'combinatorial';
[L,W,D,Di] = MeshLaplacian(X,F,lapl);
% combinatorial laplacian

q = 70;
niter = 1200;
ndisp = round( 1 + (niter-1) * linspace(0,1,q).^4 );
ndisp = unique(ndisp); q = length(ndisp);
idisp = 1;


L1 = Di*L;
% relax
tau = .8;
Xp = X;
for it=1:niter
    % smoothing displacement
    V = (L1*Xp')';
    % displacement toward sphere, grad of 1/2 * (|x|-1)^2
    U = ( sqrt(sum(Xp.^2))-1 ) .* Xp ./ sqrt(sum(Xp.^2));
    if it==1
        eta = .5 * norm(V(:)) / norm(U(:));
    end
    % grad descent
    Xp = Xp - tau*( V + eta*U  );
    % draw
    if it == ndisp(idisp)
        clf; hold on;
        plot_mesh(Xp,F,options);
        % view(60,25);
        camlight; 
        shading interp;
        % shading faceted;
        drawnow;
        saveas(gcf, [rep  'anim-' znum2str(idisp,2) '.png'] );
        idisp  = idisp+1;
    end
end


