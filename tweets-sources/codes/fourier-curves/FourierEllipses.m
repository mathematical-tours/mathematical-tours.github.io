%%
% display Fourier synthesis using sums of ellipses

name = 'triskel';
name = 'cat';

addpath('../toolbox/');
rep = MkResRep(['ellipses/' name]);


q = 512; 
I = load_image(name, q);
I = rescale(sum(I,3));

t = linspace(0,1,q);
C = contour(t,t,I,[.5,.5]);
m = C(2,1);  c = C(:,2:m+1); c = c(1,:)'+1i*c(2,:)';

% resample
p = 1024;
curvabs = @(gamma)[0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
resample1 = @(gamma,d)interp1(d/d(end),gamma,(0:p-1)'/p, 'linear');
resample = @(gamma)resample1( [gamma;gamma(1)], curvabs( [gamma;gamma(1)] ) );
c = resample(c);

% Fourier coefs
P = 50; % # frame
rmax = P;
cf = fft(c); 

t = linspace(0,1,256);
ellip = @(m,a,b,t)a*exp(2i*pi*t) + b*exp(-2i*pi*t);
Ellip = @(m,a,b)ellip(m,a,b,t);

% #coefs for the synthesis 
K = 10;

% approximated curve
f = [0:p/2, -p/2+1:-1]';
c1 = ifft( cf .* (abs(f)<=K) );

% #time
T = 100;
for i=1:T
    t = (i-1)/(T-1);
    m = cf(1)/p; % initial position
    clf; hold on;
    axis equal;  axis ij; axis off;
    % plot partial curve
    for k=1:K
        % draw an ellipse
        a = cf(1+k)/p; b = cf(end-k+1)/p;
        plot(m, 'k.' , 'MarkerSize', 20);
        plot( m + Ellip(m,a,b), 'Color', [1 1 1]*.5, 'LineWidth', 1 );
        % shift position
        m = m + ellip(m,a,b,k*t);        
    end
    plot(c1([1:end 1]), 'LineWidth',1, 'Color', [1 .3 .3]); 
    plot(c1(1:round(t*end)), 'LineWidth',2, 'Color', [1 0 0]); 
    plot(c1(max(1,round(t*end))), '.', 'MarkerSize',30, 'Color', [1 0 0]); 
    axis([0 1 0 1]); drawnow;
    saveas(gcf, [rep name '-' num2str(K) '-' znum2str(i,3) '.png']);
end

% AutoCrop(rep, [name '-'])