%%
% display Julia sets

if not(exist('cnt'))
    cnt = 1;
end

addpath('../toolbox/');
rep = MkResRep(num2str(cnt));

% center for the zoom




c = -0.4+0.6i;
z_tgt = -0.198629591421+0.144262869667i;


c = -0.624 + 0.435i;
z_tgt = 0.044016561386-0.232510349130i;

c = -0.7269+0.1889i;
z_tgt = 0.013117482253-0.184680146198i;

c = -0.835-0.2321i;
z_tgt = -0.036787892747-0.121310465626i; 

% zoom factor
rho = 1.2;

n = 256*2; % size image
niter_frac = 1000;

explore_mode = 0;

a = 1.8;




ndisp = 100;
z = 0;
for i=1:ndisp
    [Q,tx,ty] = Julia(c,z_tgt,a,n,niter_frac);
    a = a/rho;
    % z = (1-1/rho)*z + 1/rho * z_tgt;    
    %
    A = log(Q);
    clf;
    imagesc(ty,tx,A);
    colormap parula(256);
    axis image; drawnow; axis ij;
    if explore_mode && mod(i,10)==1
        [ax,ay,button] = ginput(1);
        z_tgt = ay + 1i*ax;
    end
    %
    imwrite(rescale(A)*255,parula(256),[rep 'julia-' znum2str(i,3) '.png'],'png');
end

cnt = cnt+1;