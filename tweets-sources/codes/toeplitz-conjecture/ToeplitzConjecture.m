
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);



% resample a curve using a fixed number of points
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );


name = 'elephant';
name = 'triskel';
name = {'cat' 'france'};
name = {'bunny' 'elephant'};
name = {'bunny' 'cat'};
n = 512;
p = 10000;
for i=1:length(name)
    f{i} = load_image(name{i}, n);
    f{i} = sum(f{i},3);
    f{i} = double( rescale(f{i},0,1)<.5 );
    
    
    [C,h] = contour(f{i},[.5,.5]);
    m = C(2,1);  
    c = C(:,2:m+1); c = c(1,:)'+1i*c(2,:)';
    c = (c-1)/(n-1);
    clist{i} = resample(c,p);
end


q = 100;

ntries = 5000;
tol = .1;

for it=1:q
    s = (it-1)/(q-1);
    F = (1-s)*f{1} + s*f{2};
    c = (1-s)*clist{1} + s*clist{2};    
    
    % extract square
    d = Inf; S = [];
    draw_on = 0;
    for k=1:ntries
        i = 1+floor(rand*(p));
        j = 1+floor(rand*(p));
        epsi = sign(randn);
        %
        r = abs(c(i)-c(j));
        t = angle(c(i)-c(j));
        a = c(j) + sqrt(2)*r*exp(1i*(t + epsi*pi/4) );
        b = c(j) + r*exp(1i*(t + epsi*pi/2) );
        d1 = min(abs(c-a)) + min(abs(c-b));
        if d1<d && abs(c(i)-c(j))>tol
            S = [c(j), c(i), a, b];
            d = d1;
        end
        
        if draw_on
            clf; hold on;
            plot(c, 'k', 'LineWidth', 2);
            plot(S([1:end 1]), 'r.-', 'MarkerSize', 25,'LineWidth', 2);
            axis equal; hold on;
            drawnow;
        end
        
    end

    clf; hold on;
    plot(c, 'k', 'LineWidth', 2);
    plot(S([1:end 1]), 'r.-', 'MarkerSize', 25,'LineWidth', 2);
    axis equal; hold on; axis off;
    axis ij;
    drawnow;
    mysaveas(it);
end
