%%
% Generate pixelization.

name = 'eiffel-tower-sunrise';
rep = 'results/pixelize/';
if not(exist(rep))
    mkdir(rep);
end

addpath('./images/');
base_path = '/Users/gabrielpeyre/Dropbox/github/numerical-tours/matlab/';
addpath([base_path 'toolbox_general/']);
addpath([base_path 'toolbox_signal/']);

% target size
N = 512*2;

f = load_image(name);
r = size(f,1)/size(f,2);
P = round(N/r);
f = image_resize(f, [N P 3]);
f = rescale(f);
N = size(f);

Nbr = 10000;
m_min = round( 5*(N/512) );
m_max = round( 100*(N/512) );
g = f;
remap = @(t)t(2)^2;
remap = @(t)1.5 * ( (t(1)-1/2)^2 + (t(2)-1/2)^2 );
for i=1:Nbr
    k = floor( [rand*N(1) rand*N(2)] ) + 1;
    s = ceil(m_min + rand*(m_max-m_min))+1;
    t = [(k(1)-1)/(N(1)-1), (k(2)-1)/(N(2)-1)];
    v = .5+.5*rand; v = 1;
    s = ceil(v*m_max*remap(t));
    selx = k(1):k(1)+s-1; 
    sely = k(2):k(2)+s-1; 
    if 0
    selx = mod(selx-1,N(1))+1;
    sely = mod(sely-1,N(2))+1;
    else
       selx = min(selx,N(1));
       sely = min(sely,N(2));
    end
    for q=1:3
        g(selx,sely,q) = mean(mean(f(selx,sely,q)));
    end
end

imageplot(g);
imwrite(rescale(g), [rep name '-' num2str(N) '.png'], 'png');
