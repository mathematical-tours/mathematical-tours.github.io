function F = perform_ar_synthesis(H, options)

% perform_ar_synthesis - perform motion cloud synthesis
%
%   F = perform_ar_synthesis(H, options);
%
%   Copyright (c) 2013 Gabriel Peyre

n = size(H,1);

synth2d = @(h)real(ifft2(fft2(randn(n,n)).*h));
extend = @(f)[f f(:,1); f(1,:) f(1)];


% movement
scale = getoptions(options, 'scale', 1); 
rotation = getoptions(options, 'rotation', 0); 
translation = getoptions(options, 'translation', [0 0]);
center = getoptions(options, 'center', [n/2 n/2]);
sigmat = getoptions(options, 'sigmat', 40);
p = getoptions(options, 'nbframes', 64);
pdisp = getoptions(options, 'nbframes_disp', 2*p);
%
a = fit_ar2(sigmat);


x = 1:n;
[Y,X] = meshgrid(x,x);
% moved grid
X1 = center(1) + scale*( (X-center(1))*cos(rotation) - (Y-center(2))*sin(rotation) ) + translation(1);
Y1 = center(2) + scale*( (X-center(1))*sin(rotation) + (Y-center(2))*cos(rotation) ) + translation(2);
X1 = mod(X1-1,n)+1;
Y1 = mod(Y1-1,n)+1;
% apply movement
move = @(f)interp2(1:n+1,1:n+1,extend(f),Y1,X1);

% Stop button

clf;
% Axes
ax = axes(...
    'Units','Normalized',...
    'OuterPosition', [0 0.2 1 0.8]);


f0 = zeros(n);
f1 = zeros(n);
Contrast = 50;
s = 0;
F = [];
EarlyStop = 0;
i = 0; 
randn('state', 123);
global run;
run = 1;
while run
    i = i+1;
    % rotate the noise
    W = synth2d(H);
    X1 = center(1) + ( (X-center(1))*cos(rotation*(i-1)) - (Y-center(2))*sin(rotation*(i-1)) );
    Y1 = center(2) + ( (X-center(1))*sin(rotation*(i-1)) + (Y-center(2))*cos(rotation*(i-1)) );
    X1 = mod(X1-1,n)+1;  Y1 = mod(Y1-1,n)+1;
    W = interp2(1:n+1,1:n+1,extend(W),Y1,X1);
    %
    f = W + a(1)*f0 + a(2)*f1;
    f1 = f0; f0 = f;
    s = max(s,std(f(:)));
    % scale    
    f0 = move(f0);
    f1 = move(f1);
    % display
    image( 127 + Contrast*f/s  ); 
    colormap gray(256);
    axis image; axis off;
    
    if i==p
uicontrol(...
    'Style','pushbutton', 'String', 'Stop',...
    'Units','Normalized', 'Position', [0.4 0.1 0.2 0.1],...
    'Callback', @MyCallback);
    end
    if i>=p
        set(ax, 'ButtonDownFcn', 'get(ax, ''CurrentPoint'')');
    end
    
    drawnow;
    % save
    if i<p
        F(:,:,end+1) = f;
    end
end


end

function MyCallback(a,b,c)
global run;
run = 0;
end

