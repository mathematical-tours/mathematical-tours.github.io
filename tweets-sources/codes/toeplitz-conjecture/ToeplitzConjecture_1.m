%%
% https://en.wikipedia.org/wiki/Inscribed_square_problem


addpath('../toolbox/');
addpath('../subdivision-curves/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep  'anim-' znum2str(it,3) '.png']);

%% Generate square

square = @(c,r,t)exp(2i*pi*((0:3)'/4+t))*r + c;


h4pt = @(w)[-w, 0, 1/2+w, 1, 1/2+w, 0, -w];
h = h4pt(1/16);

subdivide = @(f,h)cconvol( upsampling(f), h);
mysub = @(f)f;
for i=1:5
    mysub = @(f)mysub(subdivide(f,h));
end

k = 3;
x0 = rand(k,1) + 1i*rand(k,1);


c = [.3+.3i, .7+.5i];
c = [.5+.5i];

I = randperm(k+length(c)*4);

z = {};
q = 70;
for it=1:q
    s = (it-1)/q;       
    %
    x = x0;
    for j=1:length(c)
        delta = .1*(j-1)/(length(c)-1);
        if length(c)==1
            delta=0;
        end
        r = .05+.3*sin(pi*(s+delta)).^2;
        z{j} = square(c(j), r, s + delta );
        x = [x;z{j}];
    end
    x1 = mysub(x(I));
    %
    clf; hold on;
    plot(x1([1:end 1]), 'b-', 'LineWidth', 2);
    for j=1:length(c)
        plot(z{j}([1:end 1]), 'r.-', 'LineWidth', 2, 'MarkerSize', 30);
    end
    axis equal; axis tight; axis off;
    axis([0 1 0 1])
    drawnow;
    mysaveas(it);
end


