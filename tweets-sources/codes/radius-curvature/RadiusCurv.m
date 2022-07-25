addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


% resample a curve using a fixed number of points
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );
plotcol = @(x,y,col)surface([x(:)';x(:)'],[y(:)';y(:)'],zeros(2,length(x(:))),[col(:)';col(:)'],'facecol','no','edgecol','interp','linew',2);
closecurv = @(f) f([1:end 1]);   
    
subdivide = @(f,h)cconvol( upsampling(f), h);
h4pt = @(w)[-w, 0, 1/2+w, 1, 1/2+w, 0, -w];
h = h4pt(1/16);

diffper = @(f)( f([2:end 1])-f([end 1:end-1]) )*length(f)*2;
dotp = @(f,g)real(f.*conj(g));
normalize = @(x)x./abs(x);
        
% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 6; x = rand(k,1)+1i*rand(k,1);
eta = .02; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
q = 120;
for it=1:q
    f = x;
    for isub=1:5
      	f = subdivide(f,h);
    end        
    p = 512;
    f = resample(f,p);
    % normals
    T = normalize( diffper(f));
    N = 1i*T;
    Nkappa = diffper(T);
    kappa = dotp(Nkappa,N);
    
    % 
    clf; hold on;
    if 0
        plot(f([1:end 1]), 'k-', 'LineWidth', 2);
    else
        plotcol = @(x,y,col)surface([x(:)';x(:)'],[y(:)';y(:)'],zeros(2,length(x(:))),[col(:)';col(:)'],'facecol','no','edgecol','interp','linew',2);
        col = clamp(kappa, -u,u);
        plotcol(real(closecurv(f)), imag(closecurv(f)), closecurv(col));
        caxis([-1,1]*200);
        colormap jet(256);
    end
    % draw normals
    m = 10;
    I = round(linspace(0,p,m+1)'); I(1)=[];
    for k=1:m
        s = (k-1)/(m-1);
        i = I(k);        
        tau = .1;
        plot(f(i), '.', 'color', [s 0 1-s], 'MarkerSize', 15);
        plot([f(i) f(i)+tau*N(i)], '-', 'color', [s 0 1-s], 'LineWidth', 2);
    end
    e = .2;
    axis equal;  axis([-e,1+e,-e,1+e]); 
    set(gca, 'XTick', [], 'YTick', []); box on;
    drawnow;
    % DO HERE STUFF
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end