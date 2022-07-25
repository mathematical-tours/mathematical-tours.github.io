addpath('../toolbox/');
addpath('./mexEMD/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


stdize = @(x)(x-mean(x))/std(x);

% resample a curve using a fixed number of points
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );
closecurv = @(f) f([1:end 1]);   
%    
subdivide = @(f,h)cconvol( upsampling(f), h);
h4pt = @(w)[-w, 0, 1/2+w, 1, 1/2+w, 0, -w];
h = h4pt(1/16);

k = 4;


x0 = rand(k,1)+1i*rand(k,1);
y0 = rand(k,1)+1i*rand(k,1);
for isub=1:5
      	x0 = subdivide(x0,h);
      	y0 = subdivide(y0,h);
end        
n = 128*2;
m = n;
x0 = resample(x0,n);
y0 = resample(y0,m);

% Classical OT.
a = ones(n,1)/n;
b = ones(m,1)/m;
C = abs(x0-transpose(y0)).^2;
[cost,gamma] = mexEMD(a,b,C);
[I,J,gammaij] = find(gamma);

clf; plot_plan(x0,y0, gamma);

%%
% Now GW.

q = 20;

for it=1:q
    t = (it-1)/(q-1);
    x = x0;
    y = (1-t)*x0+t*y0;
    
    CX = abs(x-transpose(x));
    CY = abs(y-transpose(y));
    
    gamma0 = ComputeLB(CX,CY);
    if 1
        [gamma_lb,GWlb] = ComputeGW(CX,CY,gamma0);
    else
        gamma_lb = gamma0;
    end
    clf; plot_plan(stdize(x),stdize(y)*3, gamma_lb);
    axis tight;
    drawnow;

end

return;

t = .5;
x = x0;
y = (1-t)*x0+t*y0;

CX = abs(x-transpose(x));
CY = abs(y-transpose(y));

I = randperm(n);
gamma = sparse(1:n,I,ones(n,1)/n, n,n);

clf; plot_plan(stdize(x),stdize(y)*3, gamma);
axis tight;
drawnow;

clf;



return;


% init using random permutations
GWbest = +inf;
ntries = 200;
for k=1:ntries
    progressbar(k,ntries);
    I = randperm(n);
    gamma = sparse(1:n,I,ones(n,1)/n, n,n);
    [gamma,GWlist] = ComputeGW(CX,CY,gamma);
    if GWlist(end)<GWbest
        GWbest = GWlist(end);
        gamma_best = gamma;
    end
end
fprintf( 'Random/LB: %f\n', GWbest/GWlb(end) );

clf; plot_plan(stdize(x),stdize(y)*3, gamma_best);
axis tight;

