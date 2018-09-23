%%
% Test for 3D isosurface interpolation.


n = 120; 

t = linspace(-1,1,n);
[x,y,z] = ndgrid(t,t,t);

% colors
col = [1 0 0];

if 0
    F{0} = ImplicitEquation(x,y,z,'ball');
    % F0 = ImplicitEquation(x,y,z,'trumpet');
    F{0} = ImplicitEquation(x,y,z,'torus2');
    
    F{1} =    ImplicitEquation(x,y,z,'torus') .* ...
        ImplicitEquation(y,z,x,'torus') .* ...
        ImplicitEquation(z,x,y,'torus');
    F{1} = ImplicitEquation(x,y,z,'pball');
    F{1} = ImplicitEquation(x,y,z,'torus');
else
    names = {'duck' 'moomoo'};
    names = {'moomoo' 'elephant'};
    names = {'duck' 'elephant'};
    names = {'bunny' 'elephant'};
    F = {};
    for k=1:2
        [V,Fc] = read_off([names{k} '.off']);
        if strcmp(names{k}, 'elephant')
            V(1,:) = -V(1,:);
            Fc = Fc([2 1 3],:);
        end
        if strcmp(names{k}, 'bunny')
            V(3,:) = -V(3,:);
            Fc = Fc([2 1 3],:);
        end
        V = V - repmat(mean(V,2), [1 size(V,2)]);
        V = V / max(abs(V(:))) * .9;
        % signed distance function
        A = inpolyhedron(Fc',V',t,t,t);
        D = pointCloud2TDF(V,y,x,z);
        F{k} = -(2*A-1) .* D;
        % smooth a bit
        F{k} = smooth3d(F{k},2);
    end    
end

rep = MkResRep([names{1} '-' names{2}]);


q = 50;
for it=1:q
    s = (it-1)/(q-1);
    Fs = (1-s)*F{1} + s*F{2};
    col = [s 0 1-s];
    % display -- classical 
    if 1
    clf;
    p = patch( isosurface( t,t,t, Fs, 0 ) );
    isonormals( t,t,t,Fs,p );
    set(p, 'FaceColor', col, 'EdgeColor', 'none');
    box on; axis off; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    axis equal; axis([-1 1 -1 1 -1 1]);    
    view(170,-40); % view(-150,-40);
    zoom(1.1); lighting gouraud;
    camlight; drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    else
    % display -- sliced
    slide_dim = 2; slice_pos = .5;
    H = -Fs/max(Fs(:));
%    H = rescale(H);
    clf;
    VolumetricSlicing(permute(H, [3 2 1]),slide_dim, col,slice_pos);
    view(3); camlight;
    saveas(gcf, [rep 'sliced-' znum2str(it,2) '.png']);
    end
end


% AutoCrop(rep, ['sliced-']); 