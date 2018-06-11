%% Fast Marching in 2D
% This tour explores the use of Fast Marching methods in 2-D.

name = 'mountain';
name = 'bump';
name = 'piececonst';
name = 'cavern';
name = 'maze1';
name = 'obstacles';

addpath('../toolbox/');
rep = MkResRep(name);

clear options;

n = 300;
nmax = n*n;
switch name
    case 'bump'
        if not(exist('f'))
            t = linspace(0,1,n);
            [U,V] = meshgrid(t,t);
            s = .08;
            G = @(u,v)exp(-((U-u).^2+(V-v).^2)/(2*s^2));
            %
            f = zeros(n);
            while true
                axis([0 1 0 1]);
                imagesc(t,t,f');
                [a,b,button] = ginput(1);
                if button==3
                    break;
                end
                f = f+G(b,a);
            end
            f = rescale(-f, .01,1);
        end
    case 'piececonst'
        f(1:end/2,1:end/2) = 1;
        f(end/2+1:end,1:end/2) = .75;
        f(1:end/2,end/2+1:end) = .25;
        f(end/2+1:end,end/2+1:end) = .5;
    otherwise
        f = rescale( sum(load_image(name, n),3) );
        if strcmp(name, 'maze1')
            f = double(f>.5);
        elseif  strcmp(name, 'obstacles')
            f(f<.15) = 0;
        end
        nmax = sum(f(:)>0);
        f = rescale(f,1e-8,1);
end

% propagates slowly in area of large values


%%
% Display the image.

clf;
imageplot(f);

%%
% Define start and end points \(x_0\) and \(x_1\) (note that you can use your own points).

if not(exist('X'))
X = [];
% click selection
it = 0;
clf; hold on;
imagesc(1:n,1:n,f'); axis equal;
while true
    axis([1 n 1 n]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    X(:,end+1) = [a;b];
end
end

options.nb_iter_max = Inf;
options.end_points = [];
if 0
[Dfull,S] = perform_fast_marching(f, X, options);
%
clf;
hold on;
imagesc(1:n,1:n, convert_distance_color(Dfull,f,256, max(Dfull(:))));
plot(X(2,:), X(1,:), 'b.', 'MarkerSize', 25); 
end

% display progressive propagation.
q = 50; % #display
niter = round(linspace(0,1,q)*nmax);
for i=1:q
    options.nb_iter_max = niter(i);
    options.end_points = [];
    [D,S] = perform_fast_marching(f, X, options);  
    
    if 0
    clf; hold on;
    imageplot( convert_distance_color(D,f) );
    plot(X(2,:), X(1,:), 'b.', 'MarkerSize', 25);
    end
    
    % display quantized colormap
    D(S==0) = Inf;
    D1 = D; % D1(D1==Inf) = D1(D1~=Inf); % NaN
    vmax = max(D1(:));
    r = 25; % #levellines
    clf; hold on;
    [Dc,vmax] = convert_distance_color(D,f,r);
    imagesc(1:n,1:n, Dc );
    vlist = linspace(0,vmax,r); vlist(end)=vlist(end)-1e-4;
    contour(1:n,1:n,D1,vlist, 'k');
    axis image; axis off;

    saveas(gcf, [rep name '-fm-' znum2str(i,3) '.png' ]);
end


% AutoCrop(rep, [name '-fm-'])


