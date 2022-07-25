addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

k = 256/2;
n = 256*2;

name = 'shannon';
name = 'phantom';
f0 = rescale( load_image(name, 2*k) );

% f0 = f0.*(t*t');

sel0 = n/2-k+1:n/2+k;
f = zeros(n);
f(sel0,sel0) = f0;

F = abs(fft2(f));

g = real( ifft2(F) );
imageplot(clamp(g));

ProjPhase = @(G,F)real(ifft2( G./abs(G) .* F  ));
ProjPhase = @(g,F)ProjPhase(fft2(g),F);
%
rProjPhase = @(F,f)2*ProjPhase(F,f)-F;
rProjZeros = @(F,sel)2*ProjZeros(F,sel)-F;

meth = 'pocs';
meth = 'dr1';
meth = 'dr';

% imwrite(clamp(f(sel0,sel0)), [rep  'input.png']);

sel = find(f(:)<1e-3);

niter = 300; % for DR
mu = 1; % for DR
g = rand(n);
g1 = rand(n);
Err = [];
for it=1:niter
    switch meth
        case 'pocs'
            g = ProjPhase(g,F);
            g = ProjZeros(g,sel);        
        case 'dr'
            g1 = (1-mu/2)*g1 + mu/2*rProjPhase( rProjZeros(g1,sel), F);
            g = ProjZeros(g1,sel);
        case 'dr1'
            g1 = (1-mu/2)*g1 + mu/2*rProjZeros(rProjPhase(g1,F),sel);
            g = ProjPhase(g1,F);
    end
    Err(end+1) = norm(abs(fft2(g))-F, 'fro');
    imwrite(clamp(g(sel0,sel0)), [rep  'anim-' znum2str(it,3) '.png']);
    if 1 % mod(it,2)==1
        clf; imageplot(clamp(g(sel0,sel0)));
        drawnow;
    end
end