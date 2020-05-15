%%
% Display the core of a polygon

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


% click selection
x = {};
clf; hold on;
col = 'k';
for num=1:2
    x{num} = [];
    while true
        if num==2 && length(x{1})==length(x{2})
            break;
        end
        axis equal; axis([0 1 0 1]);
        box on; set(gca, 'Xtick', [], 'Ytick', []);
        [a,b,button] = ginput(1);
        if button==3
            break;
        end        
        plot(a,b, '.', 'color', col, 'MarkerSize', 15);
        x{num}(end+1) = a + 1i*b;
    end
    col = 'r';
    x{num} = x{num}(:);
end

m = length(x{1});

n = 512*2;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;
dotp = @(a,b)real(a.*conj(b));


% anim
q = 50;
for it=1:q
    s = (it-1)/(q-1);
    z = (1-s)*x{1}+s*x{2};
    col = [s 0 1-s];
    clf; hold on;
    M = zeros(n); C = ones(n);
    for i=1:m
        a = z(i); b = z(mod(i,m)+1);
        u = 10;
        I = dotp(Z-a,1i*(b-a))>=0;
        M(I) = M(I)+.05;   
        C = C .* I;
    end
    % C = (M==max(M(:)));
    CA = ones(n,n,3);
    for j=1:3
        u = ones(n); u(C==1) = col(j);
        CA(:,:,j) = u;
    end
    imagesc(t,t,1-M', 'AlphaData', .6);
    colormap gray(256); caxis([0 1]); 
    fill(real(z), imag(z), .2*col+.8*[1 1 1], 'EdgeColor', col, 'FaceAlpha', .5); 
    % alpha(.2);    
    for i=1:m
        a = z(i); b = z(mod(i,m)+1);
        u = 10;
        plot( [a*u + b*(1-u),-a*u + b*(1+u)], 'color', [1 1 1]*.6 );
    end   
    plot(z([1:end 1]), 'color', col,'LineWidth', 2); 
    h = imagesc(t,t,permute(CA, [2 1 3])); set(h,'AlphaData',C'==1);
    axis equal; axis([0 1 0 1]); box on;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    mysaveas(it);
end