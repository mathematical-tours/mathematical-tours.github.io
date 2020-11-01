%%
% Wave vs. Heat vs. Schrodinger equation in the Fourier domain


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% step size
n = 1024;

% Laplacian in frequency domain
omega = 2*[0:n/2, -n/2+1:-1]'/n;
OmSq = omega.^2; 

% initialization

x = linspace(-1,1,n)';
s = .1;

gauss = @(m,s)exp( -(x-m).^2/(2*s^2) );
f0 = abs(x)<=.4;
f0 = gauss(0,s);
f0 = gauss(-.3,.08) - .7*gauss(.3,.05);
f0 = gauss(-.3,.1) - .7*gauss(.3,.08);

[Y,X] = meshgrid(t,t);
F0 = fft(f0);

name = 'biheat';
name = 'biwave';
name = 'wave';
name = 'schrodinger';
name = 'heat';

mydisp = @(f)real(f);
switch name
    case 'heat'
        % heat, df/dt = Delta(f)
        Tmax = 5*150^2;
    case 'biheat'
        % df/dt = Delta^2(f)
        Tmax = 100^4;
    case 'wave'
        % heat, df^2/dt = Delta(f)
        Tmax = 500;     
    case 'biwave'
        % heat, df^2/dt = Delta^2(f)
        Tmax = 10^2;  
    case 'schrodinger'
        % heat, i*df^2/dt = Delta(f)
        Tmax = 50000; 
        mydisp = @(f).05 + 1.95*abs(f).^2-1;
end

q = 100; % #frames

% svg run
F = [];
for it=1:q  
    % 
    s = (it-1)/(q-1);
    t = s*Tmax;  % time
    switch name
        case 'heat'
            Ft = F0 .* exp( -OmSq*t );
        case 'biheat'
            Ft = F0 .* exp( -OmSq.^2*t );
        case 'wave'
            Ft = F0 .* exp( 2i*pi*sqrt(OmSq)*t );
        case 'schrodinger'
            Ft = F0 .* exp( 2i*pi*(OmSq)*t );
        case 'biwave'
            Ft = F0 .* exp( 2i*pi*OmSq*t );
    end
    F(:,it) = ifft(Ft);
    
end


% svg run
for it=1:q  
    % 
    s = (it-1)/(q-1);
    clf; hold on;
    cnt = 0;
    for jt=max(1,it-20):it
        cnt = cnt+1;
        r = (cnt-1)/20;
        plot(mydisp(F(:,jt)), 'color', (1-.5*r) + .5*r*[s 0 1-s], 'LineWidth', 2);
    end
    plot(mydisp(F(:,it)), 'color', [s 0 1-s], 'LineWidth', 3);    
    axis([0,n,-1,1 ]);
    box on;
    set(gca, 'XTick', [], 'YTick', []);    
    drawnow;
    %
    mysaveas(it);
end


% AutoCrop(rep, [name '-']);
% convert heat-*.png heat.gif
% convert wave-*.png wave.gif
% convert biheat-*.png biheat.gif
% convert biwave-*.png biwave.gif