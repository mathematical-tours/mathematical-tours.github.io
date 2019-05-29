function [g,f] = pixelize(f, height, options)

options.null = 0;
Nbr = getoptions(options, 'nsamples', 10000);


remap = getoptions(options, 'remap', @(t)t(2)^2);
rep = getoptions(options, 'rep', []);

r = size(f,1)/size(f,2);
P = round(height/r);
f = image_resize(f, [height P 3]);
f = rescale(f);
N = size(f);

q = 50; % #display
ndisp = round(linspace(1,Nbr,q));
kdisp = 1;

m_max = round( 100*(height/512) );
g = f*0+1;
for i=1:Nbr
    k = floor( [rand*N(1) rand*N(2)] ) + 1;
    t = [(k(1)-1)/(N(1)-1), (k(2)-1)/(N(2)-1)];
    s = ceil(m_max*remap(t));
    selx = k(1):k(1)+s-1;
    sely = k(2):k(2)+s-1;
    selx = min(selx,N(1));
    sely = min(sely,N(2));
    for q=1:3
        g(selx,sely,q) = mean(mean(f(selx,sely,q)));
    end
    if ndisp(kdisp)==i
        clf;
        imageplot(g); drawnow;
        if not(isempty(rep))
            imwrite(rescale(g), [rep 'anim-' znum2str(kdisp,2) '.png'], 'png');
        end
        kdisp = kdisp + 1;
    end
end


end