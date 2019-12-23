%%
% Display B and E fields

addpath('../toolbox/');
rep = MkResRep();
flat = @(x)x(:);

randn('state', 2124);
rand('state', 2124);

field = 'magnetic';
field = 'electric';

animation = 'merge';
animation = 'straigth';


n = 25;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;


m = 300;
ti = linspace(0,1,m);
[Yi,Xi] = meshgrid(ti,ti);
Zi = Xi+1i*Yi;

% random mono/dipole locations
k = 12;
k = 4;
k = 2;


% orientations
delta = rand(1,1,k) + 1i*rand(1,1,k);
delta = delta./abs(delta);

p0 = rand(1,1,k) + 1i*rand(1,1,k);
s = reshape((-1).^(1:k), [1 1 k]);

% speed 
eta = .05;
v = randn(1,1,k) + 1i*randn(1,1,k);
v = eta*v./abs(v);


switch field
    case 'electric'
        Ef = @(p,delta)sum( s .* (Z - p) ./ abs(Z-p).^3 ,3);
    case 'magnetic'
        Ef = @(p,delta)sum( ...
                3 * (Z - p) .* real( (Z-p).*conj(delta) ) ./ abs(Z-p).^5 - ...
                delta ./ abs(Z-p).^3 ,3 );
end


% rotation speed
eta_r = eta;
eta_r = 0.02;
r = exp( (-1).^(1:k) * 2i*pi*eta_r );
r = reshape(r, [1 1 k]);

p = p0;

q = 75;
for it=1:q
    % field
    E = Ef(p,delta); 

    % plot
    U = E./abs(E);
    clf; hold on;
    % imagesc(ti,ti, log(abs(PE')));
    % caxis([0 10]);    
	quiver(t,t,imag(U), real(U), 'k', 'filled', 'LineWidth', 1, 'AutoScaleFactor', .7); 
    switch field
        case 'electric'
            z = flat(p(s<0)); plot( imag(z), real(z), 'r.', 'MarkerSize', 25 );
            z = flat(p(s>0)); plot( imag(z), real(z), 'b.', 'MarkerSize', 25 );
        case 'magnetic'
        	quiver(imag(p(:)),real(p(:)),.07*imag(delta(:)),.07*real(delta(:)), 'm', 'filled', 'LineWidth', 3, 'AutoScale', 0); 
    end
    g = .07;
    axis([-g 1+g -g 1+g]); 
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
    box on;  axis off; axis ij;
    drawnow;
    
    switch animation
        case 'straigth'
            % advance
            p = p + v;
            % rotate
            delta = delta .* r;
            % reflexion on boundary
            I = find( real(p)<0 | real(p)>1 );
            v(I) = -real(v(I)) + 1i*imag(v(I));
            I = find(  imag(p)<0 | imag(p)>1 );
            v(I) = real(v(I)) - 1i*imag(v(I));
        case 'merge'
            T = (it-1)/(q-1);
            if 0
                c = ( p0(1:2:end) + p0(2:2:end) )/2;
                cd = ( p0(1:2:end) - p0(2:2:end) )/2;
                p(1:2:end) = c + (1-T) * cd;
                p(2:2:end) = c - (1-T) * cd;
            else
                p(2:2:end) = (1-T)*p0(2:2:end) + T*p0(1:2:end);
            end
    end
    
    
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end
