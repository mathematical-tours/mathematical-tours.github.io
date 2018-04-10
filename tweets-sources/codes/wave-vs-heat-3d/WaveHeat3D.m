%%
% Wave and heat in 3D.


model = 'heat'; 
model = 'wave'; 

% 'horse' 'shears' 
name = 'duck';
name = 'elephant';
name = 'moomoo';
name = 'cube';
name = 'gaussians';

addpath('../toolbox/');
rep = MkResRep(name);


n = 120; 
x = linspace(-1,1,n);
[X,Y,Z] = ndgrid(x,x,x);

pov = 3;
switch name
    case 'gaussians'
        %
        M = {[.5 .3 .8] [.5 .7 .2] [.8 .5 .5]};
        M = {[.6 .5 .5]};
        w = [1 .8 .7];
        %
        s = .08;
        f0 = zeros(n,n,n);
        for k=1:length(M)
            z = 2*M{k}-1;
            DD = (X-z(1)).^2 + (Y-z(2)).^2 + (Z-z(3)).^2;
            f0 = f0 + w(k) * exp(-DD/(2*s^2));
        end
    case 'cube'
        f0 = max( max(abs(X), abs(Y)), abs(Z) )<=.6;
    otherwise
        [V,F] = read_off([name '.off']);
        V = V - repmat(mean(V,2), [1 size(V,2)]);
        V = V / max(abs(V(:))) * .9;
        switch name
            case 'moomoo'
                pov = [0,90];
            case 'duck'
                pov = [-130,-15];
            case 'elephant'
                pov = [5,90];
        end
        f0 = inpolyhedron(F',V',x,x,x);
end



% Laplacian in frequency domain
u = 2*[0:n/2, -n/2+1:-1]'/n;
[X,Y,Z] = meshgrid(u,u,u);
OmSq = X.^2+Y.^2+Z.^2; 

switch model
    case 'heat'
        % heat, df/dt = Delta(f)
        Tmax = 10^2; % elephant
        Tmax = 30^2; % all
    case 'wave'
        % heat, df^2/dt = Delta(f)
        Tmax = 40;     
end


q = 50; % #frames
slide_dim = 0; % mean levelset
slide_dim = 2; % cut along an axis

F0 = fftn(f0);
for i=1:q  
    % 
    t = ((i-1)/(q-1) + .01 )*Tmax;
    switch model
        case 'heat'
            Ft = F0 .* exp( -OmSq*t );
        case 'wave'
            Ft = F0 .* exp( 2i*pi*sqrt(OmSq)*t );
    end
    f = real( ifftn(Ft) );
    
    % display
   	color = [(i-1)/(q-1) 0 1-(i-1)/(q-1)];
    if slide_dim==0
        clf;
        VolumetricRendering(f, color);
    else
        slice_pos = .5;
        if strcmp(name, 'moomoo')
            slice_pos = .45;% for mooho
        end
        clf;
        VolumetricSlicing(f,slide_dim, color,slice_pos);
    end
    view(pov); camlight;

    %
    saveas(gcf, [rep model '-' num2str(slide_dim) '-' znum2str(i,2) '.png'], 'png');
end
axis tight;



% AutoCrop(rep, [model '-' num2str(slide_dim) '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif
