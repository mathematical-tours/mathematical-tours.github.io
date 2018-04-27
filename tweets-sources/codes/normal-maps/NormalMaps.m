%%
% Render normal maps using phong model

name = 'stone';
name = 'shapes';
name = 'money';
name = 'crossnrm';
name = 'texture';
name = 'lion';
name = 'rocks';
addpath('../toolbox/');
rep = MkResRep(name);

% normalize to unit norm
normalize = @(f)f ./ repmat( sqrt(sum(f.^2,3)), [1 1 3] );
rm = @(u) repmat( reshape(u,[1 1 3]), [n n 1] );

% load the map
n = 256;
f = load_image(name, n);
f = clamp(f/128-1,-1,1);
f = normalize(f);

clf; imagesc((f+1)/2);
imwrite((f+1)/2, [rep 'original.png']);

% AutoCrop(rep, [name '-fm-'])

% light position
l = [.5 .5 7];
% View position
v = [5 3 6];

model = 'dull';
model = 'diffuse';
model = 'shiny';

col_g = [1 .5 0];
col_ph = [1 1 1];
switch model
    case 'diffuse'
        k = [1 0];
        alpha = 1;
    case 'shiny'
        k = [1 1];
        alpha = 10;
    case 'dull'
        k = [1/2 1/2];
        alpha = 4;
end

q = 50; % #frame
for i=1:q
    t = (i-1)/q;
    % position of light
    r = 10;
    l = [r*cos(2*pi*t), r*sin(2*pi*t), 5];
    %
    [R,G,Ph] = RenderNormalMap(f,l,v, k, col_g, col_ph, alpha);
    % display
    clf;
    imagesc(R);
    axis image; axis off;
    caxis([0 1]);
    colormap gray(256);
    drawnow;
    % save
    imwrite(R, [rep model '-' znum2str(i,2) '.png']);
end
