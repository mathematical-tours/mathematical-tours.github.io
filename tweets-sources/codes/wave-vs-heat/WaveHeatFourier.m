%%
% Wave vs. Heat equation in the Fourier domain


addpath('../toolbox/');
rep0 = MkResRep();
rep = [rep0 '/fft/'];
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% step size
n = 256;

% Laplacian in frequency domain
t = 2*[0:n/2, -n/2+1:-1]'/n;
[Y,X] = meshgrid(t,t);
OmSq = X.^2+Y.^2; 

% initialization

t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
%
Z = [.7+.5i; .2+.8i; .4+.2i];
w = [1 -.6 .8];
%
s = .015;
f0 = zeros(n);
for k=1:length(Z)
    z = Z(k);
    f0 = f0 + w(k) * exp(-( (X-real(z)).^2+(Y-imag(z)).^2)/(2*s^2));
end
% fft 
F0 = fft2(f0);

name = 'biheat';
name = 'heat';
name = 'wave';
name = 'biwave';

switch name
    case 'heat'
        % heat, df/dt = Delta(f)
        Tmax = 150^2;
    case 'biheat'
        % df/dt = Delta^2(f)
        Tmax = 100^4;
    case 'wave'
        % heat, df^2/dt = Delta(f)
        Tmax = 50;     
    case 'biwave'
        % heat, df^2/dt = Delta^2(f)
        Tmax = 10^2;     
end

q = 50; % #frames

% svg run
for i=1:q  
    % 
    t = (i-1)/(q-1)*Tmax;  % time
    switch name
        case 'heat'
            Ft = F0 .* exp( -OmSq*t );
        case 'biheat'
            Ft = F0 .* exp( -OmSq.^2*t );
        case 'wave'
            Ft = F0 .* exp( 2i*pi*sqrt(OmSq)*t );
        case 'biwave'
            Ft = F0 .* exp( 2i*pi*OmSq*t );
    end
    f = real( ifft2(Ft) );
    
    % display
    v = max(abs(f(:)));
    g = f/v;
    r = 16; % #levellines
    clf; hold on;
    u = linspace(0,1,n);
    imagesc(u,u,g);
    contour(u,u,g,linspace(-1,1,r), 'k');
    colormap(jet(r-1));
    caxis([-1 1]);
    axis image; axis off; axis ij;
    %
    saveas(gcf, [rep name '-' znum2str(i,2) '.png'], 'png');
end
axis tight;


% AutoCrop(rep, [name '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif
% convert biheat-*.png biheat.gif
% convert biwave-*.png biwave.gif