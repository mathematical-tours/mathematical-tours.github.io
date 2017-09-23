%%
% Test for denoising 1D function.

name = 'hibiscus';
F = imread([name '.png']);

F = sum(F,3)/255;
clf; imagesc(F); axis image; colormap gray;

% extract 1-D trace
f= F(:,end/2);
n = length(f);
clf; plot(f); axis tight;

% make it periodic 
t = linspace(0,1,n)';
f = f - f(1) + (f(1)-f(end))*t;
plot(f); axis tight;

s_list = linspace(1e-6, .1, 20);
clf; hold on;
for i=1:length(s_list)
    s = s_list(i);
    % gaussian filter
    t = [0:n/2, -n/2+1:-1]'/n;
    h = exp(-t.^2/(2*s^2)); h = h/sum(h);
    % convolve
    g = ifft(fft(f).*fft(h));
    % plot
    u = (i-1)/(length(s_list)-1);
    plot(linspace(0,1,n),g,'Color', [1-u 0 u]);
end
axis tight; axis on; box on;