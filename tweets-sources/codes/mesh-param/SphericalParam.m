
addpath('../toolbox/');
rep = MkResRep('spherical');

modep = 'spherical';
modep = 'disk';

name = 'elephant';
[X,F] = read_off([name '.off']);
n = size(X,2);
clear options; options.name = name;

% center
X = X - repmat( mean(X,2), [1 n] );

clf; hold on;
plot_mesh(Xp,F,options);
camlight;


% Laplacian with cotan weights
[L,W,D,Di] = MeshLaplacian(X,F);

switch modep
    case 'spherical'
        L1 = Di*L;
        % relax
        niter = 10000;
        tau = 1;
        % project on sphere
        Xp = X ./ repmat( sqrt(sum(X.^2,1)), [3 1] );
        for i=1:niter
            Xp = Xp - tau*(L1*Xp')';
            % project on sphere
            Xp = Xp ./ repmat( sqrt(sum(Xp.^2,1)), [3 1] );
            % draw
            if mod(i,20)==1
                clf; hold on;
                plot_mesh(Xp,F,options);
                camlight; shading faceted;
                drawnow;
            end
        end
    case 'disk'
        % cut in two the surface
        s = 1; t = -.11;
        V = X(s,:)>t; % to keep 
        % faces to keep 
        I = find( V(F(1,:)) & V(F(2,:)) & V(F(3,:)) );
        F = F(:,I);
        %
        clf; hold on;
        plot_mesh(X,F,options);
        camlight;
        drawnow;
        % boundary
        options.verb = 0;
        B = compute_boundary(F, options);
        p = length(B);
        %
        t = pi*1.05 + linspace(0,2*pi(),p+1); t(p+1) = [];
        Z = [cos(t); sin(t)];
        % fixed positions Laplacian
        L1 = L;
        L1(B,:) = 0;
        for i=1:length(B)
            L1(B(i),B(i)) = 1;
        end
        % RHS
        R = zeros(2, n);
        R(:,B) = Z;        
        % Solve linear system
        Y = (L1 \ R')';
%        Xp = [Y(1,:); zeros(1,n); Y(2,:)];
        Xp = [zeros(1,n); Y(1,:); Y(2,:)];
end

% display anim sphere -> elephant
q = 70; 
for i=1:q
    t = (i-1)/(q-1);
    X1 = (1-t)*X + t*Xp;
    options.face_vertex_color = .4 + .6 * repmat([t 0 1-t], [n 1]);
    X1(1,:) = X1(1,:) + 1e-3*randn(1,n);
    clf; hold on;
    plot_mesh(X1,F,options);
    % view(120,-8);
    view(40,0)
    shading faceted;
    camlight; 
    drawnow;
    saveas(gcf, [rep modep '-' znum2str(i,2) '.png'] );
end
% AutoCrop(rep, [modep '-']); 
