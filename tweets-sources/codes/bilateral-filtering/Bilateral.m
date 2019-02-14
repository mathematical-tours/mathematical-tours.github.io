%%
% Test for bilateral filtering.

addpath('../toolbox/');
addpath('bilateral-toolbox/');
rep = MkResRep();


name = 'hibiscus';
n = 512/2; 
f = load_image(name,n);
f = clamp(f/255);
d = size(f,3);

% sigmar        : width of range Gaussian
% sigmas        : width of spatial Gaussian
sigmar = 40;
sigmas = 3;

q = 50;


sigmas = linspace(6,6,q);
sigmar = linspace(5,100,q);

sigmas = linspace(.1,50,q);
sigmar = linspace(40,40,q);

for it=1:q
    g = mybil(f,sigmar(it), sigmas(it));
    imageplot(g); 
    drawnow;
    imwrite(g, [rep name '-' znum2str(it,2) '.png']);
end