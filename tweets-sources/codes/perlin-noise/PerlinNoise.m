%%
% Gaussian texture synthesis 

addpath('../toolbox/');
rep = MkResRep();

t = [0:n/2,-n/2+1:-1]/n;
[Y,X] = meshgrid(t,t);
Gauss = @(m,s)exp(-( (X-m(1)).^2 + (Y-m(2)).^2 ) / (2*s^2));




s = .03; % spread in Fourier
om = .5; 
r = 0;



s = .02; % spread in Fourier
om = .3; 
r = .06;


q = 40;
s_list = linspace(.02,.01,q);
om_list = linspace(pi*.3,0,q);
r_list = linspace(.06,0,q);


s_list = linspace(.03,.01,q);
om_list = linspace(pi*.1,pi*.7,q);
r_list = linspace(.04,.1,q);




U = fft2(randn(n));

for it=1:q
    t = (it-1)/(q-1);
    
    s = s_list(it); % spread in Fourier
    om = om_list(it);
    r = r_list(it);
    
    
    P = Gauss( r*[cos(om) sin(om)],s );
    f = real( ifft2( U.*P ) );
    f = (f-mean(f(:)))/std(f(:));
    
    % f = abs(f);
    a = 2.5;
    f = rescale( clamp(f,-a,a) );
    
    % draw power spectrum
    sc = 3;
    P = Gauss( sc*r*[cos(om) sin(om)],sc*s ) + Gauss( -sc*r*[cos(om) sin(om)],sc*s );
    P = fftshift(P); P = P/max(P(:));
    %
    rc = 12; % #levelsets
    m = linspace(0,1,rc-1)';
    cm = m*[t 0 1-t] + (1-m)*[1 1 1];
    
    %
    clf; hold on;
    imagesc(P);
    contour(P,linspace(0,1,rc), 'k');
    colormap(cm);
    caxis([0 1]);
    axis image;  box on;
    set(gca, 'XTick', [], 'YTick', []);
    saveas(gcf, [rep 'spectrum-' znum2str(it,2) '.png']);

    
    m = linspace(0,1,256)';
    cm = m*[t 0 1-t] + (1-m)*[1 1 1];
    g = apply_colormap(f,cm);
    clf;
    imageplot(g);
    drawnow;
    imwrite(g, [rep 'synthesis-' znum2str(it,2) '.png']);

end

% AutoCrop(rep, 'spectrum');