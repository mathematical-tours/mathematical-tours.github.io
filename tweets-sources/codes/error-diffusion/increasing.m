%%
% Error diffusion for dithering.

name = 'shannon';
name = 'shannon';
method = 'simple';
method = 'floyd';

addpath('../toolbox/');
rep = MkResRep([name '-' method]);

q = 50;
nlist = round(linspace(16,128*2,q));

for it=1:q
    
    n = nlist(it);
    f = load_image(name,n);
    f = sum(f,3);
    [~,I] = sort(f(:)); f(I) = linspace(0,1,n*n);
    
    g = f;
    for i=2:n-1
        for j=2:n-1
            e = g(i,j) - (g(i,j)>.5);
            g(i,j) = g(i,j)>.5;
            switch method
                case 'simple'
                    g(i+1,j) = g(i+1,j) + e/2;
                    g(i,j+1) = g(i+1,j) + e/2;
                otherwise
                    e = e/16;
                    g(i,j+1) = g(i,j+1)     + e*7;
                    g(i+1,j+1) = g(i+1,j+1) + e;
                    g(i+1,j)   = g(i+1,j)   + e*5;
                    g(i+1,j-1) = g(i+1,j-1) + e*3;
            end
        end
    end
    
    G = imresize(g(2:n-1,2:n-1),'nearest','OutputSize',[256 256]*2)>.5;
    clf; imageplot(G); drawnow;
    imwrite(G,[rep 'anim-' znum2str(it,2) '.png']);
end