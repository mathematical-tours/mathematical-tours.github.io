%%
% Error diffusion for dithering.

name = 'shannon';
addpath('../toolbox/');
rep = MkResRep(name);

n = 160;
f = load_image(name,n);
f = sum(f,3);
[~,I] = sort(f(:)); f(I) = linspace(0,1,n*n);

q = 50; 
ndisp = round(linspace(1,(n-1)*(n-1),q));
kdisp = 1;

g = f;
it = 1;
for i=1:n-1
    for j=1:n-1
        e = g(i,j) - (g(i,j)>.5);
        g(i,j) = g(i,j)>.5;
        g(i+1,j) = g(i+1,j) + e/2;
        g(i,j+1) = g(i+1,j) + e/2;
        if ndisp(kdisp)==it
            h = clamp(g(1:end-1,1:end-1));
            clf; imageplot(h); drawnow;
            imwrite(h,[rep 'anim-' znum2str(kdisp,2) '.png']);
            kdisp = kdisp+1;
        end
        it = it+1;
    end    
end