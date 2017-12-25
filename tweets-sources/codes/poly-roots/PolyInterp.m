%%
% Display roots of polynomials

intmode = 'linear';
intmode = 'circle';

if not(exist('test'))
    test=0;
end
test = test+1;

rep = ['../results/poly-interp/' num2str(test) '-' intmode '/'];
[~,~] = mkdir(rep);

d = 7; % degree


% click and play
for k=1:2
    clf; hold on; 
    if k>1
        plot(z0{1}, 'ko');
    end
    z0{k} = [];
    for i=1:d
        axis([-1 1 -1 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if button==3
            break;
        end
        z0{k}(end+1) = a+1i*b;
    end
    z0{k} = z0{k}(:);
    % coeffs from the roots    
    x = sym('x');
    P0{k} = double( coeffs(expand(prod(x-z0{k}))) );
    P0{k} = P0{k}(end:-1:1);
    P0{k} = P0{k}/P0{k}(end);
end


% grid for vizualization
q = 201;
u = 1.1;
x = linspace(-u,u,q);
[Y,X] = meshgrid(x,x);
Z = X+1i*Y;



m = 50; % #anim
nls = 9;
for i=1:m
    % interp coeffs
    t = (i-1)/(m-1);
    switch intmode
        case 'circle'
            theta = ( 1 + exp(2i*pi*t) )/2;
        case 'linear'
            theta = t;
    end
    p = P0{1}*theta + P0{2}*(1-theta);
    
    % p = p/p(end);
    
    % display interpolating curve
    t = linspace(0,1,100);
    clf; hold on;
    for s=1:d
        switch intmode
            case 'circle'
                theta = ( 1 - exp(2i*pi*t) )/2;
            case 'linear'
                theta = t;
        end
        U = P0{1}(s)*theta + P0{2}(s)*(1-theta);
        plot( real(U), imag(U), 'LineWidth', 2, 'color', [(s-1)/(d-1) 0 1-(s-1)/(d-1)] );
    end
    plot(p, 'k.', 'MarkerSize', 30);
    axis tight; axis equal; box on; 
    set(gca, 'FontSize', 20); % set(gca,'XTick',[-1 0 1], 'YTick',[-1 0 1]);
    saveas(gcf, [rep 'roots-' num2str(i) '.png'], 'png');
    
    % display root   
    F = polyval(p,Z);
    R = roots(p);
    pvmax = 2;
    %
    clf; hold on;
    imagesc(x,x,abs(F));
    contour(x,x,(abs(F)), linspace(0,pvmax,nls), 'k', 'LineWidth', 2);
    % remove those outside
    I = find( abs(imag(R))<u & abs(real(R))<u  );
    plot(imag(R(I)),real(R(I)), 'r.', 'MarkerSize', 20);
    colormap(parula(nls-1)); caxis([0 pvmax]);
    axis([-u u -u u]);
    axis off; axis equal;
    drawnow;
    saveas(gcf, [rep 'anim-' num2str(i) '.png'], 'png');
end
