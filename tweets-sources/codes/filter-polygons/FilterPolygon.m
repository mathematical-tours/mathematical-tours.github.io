
name = 'rotating';
name = 'rotatcv';
name = 'rotatdv';
name = 'lowpass';
name = 'shift';

addpath('../toolbox/');
rep = MkResRep(name);

% click and play
if not(exist('y'))
    y0 = [];
    clf; hold on;
    while true
        axis equal; axis([-1 1 -1 1]);
        box on; set(gca, 'Xtick', [], 'Ytick', []);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 20);
        if button==3
            break;
        end
        y0(:,end+1) = a + 1i*b;
    end
end
% y0 = y0-mean(y0);
n = length(y0);
y0 = y0(:);

% circular conv
cconv = @(x,y)ifft(fft(x).*fft(y));

delta = zeros(n,1); delta(1) = 1;
t = [0:ceil(n/2), -floor(n/2)+1:-1]';

switch name
    case 'lowpass'
        h0 = [2 1 1]; h0 = h0/sum(h0);
        h = zeros(n,1); k = floor(length(h0)/2);
        h([1:k+1, end-k+1:end]) = h0;
        % slow down
        sl = .3;
    case 'shift'
        h0 = [0 1 1]; h0 = h0/sum(h0);
        h = zeros(n,1); k = floor(length(h0)/2);
        h([1:k+1, end-k+1:end]) = h0;
        % slow down
        sl = 0;
    case 'rotatcv'
        H = 1./(abs(t)+1) .* exp(2i*pi*rand(n,1)); 
        H(1) = exp(2i*pi/3);
%        H([2 end]) = abs(H([2 end])) * ;
        h = ifft(H);
        % slow down
        sl = .9;
    case 'rotatdv'
        rand('state', 123);
        H = 3./(abs(t)+1) .* exp(2i*pi*rand(n,1)); 
        H(1) = exp(2i*pi/3);
%        H([2 end]) = abs(H([2 end])) * ;
        h = ifft(H);
        % slow down
        sl = .85;
        sl = .96;
    case 'rotating'
        H = 1./(abs(t)) .* exp(2i*pi*rand(n,1)); 
        H(1) = exp(2i*pi/3);
        H  = H./abs(H);
%        H([2 end]) = abs(H([2 end])) * ;
        h = ifft(H);
        % slow down
        sl = 0;
end

% plot filter

h = sl*delta + (1-sl)*h;

q = 60;
y = y0; Y = [];
for it=1:q
    s = (it-1)/(q-1);
    % previous
    clf; hold on;
    for j=1:it-1
        r = (j-1)/(q-1);
        plot(Y([1:end 1],j), 'color', [r 0 1-r]*.3 + [1 1 1]*.7, 'LineWidth', 2);
    end
    plot(y([1:end 1]), 'color', [s 0 1-s], 'LineWidth', 3);
    axis equal; axis([-1 1 -1 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    drawnow;
    Y(:,end+1) = y;
    y = cconv(y,h);
    % saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'anim-') 