%%
% Test for min-max games optimization.

addpath('../toolbox/');
rep = '../results/min-max/';
[~,~] = mkdir(rep);

q = 501;
t = linspace(-1,1,q);
[Y,X] = meshgrid(t,t);

a = 2;
F = a*X.*Y+X.^2-Y.^2;
GradF = @(z)( a*imag(z)+2*real(z) ) - 1i * ( a*real(z)-2*imag(z) );


if not(exist('Z'))
%
clf; hold on;
% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(t,t,F');
contour(t,t,F',linspace(min(F(:)), max(F(:)),r), 'k');
colormap(parula(r-1));
caxis([min(F(:)) max(F(:))]);
axis image; axis off;
Z = [];
while true
    axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    Z(end+1) = a+1i*b;
end
end

m = 10;
Z = exp( 2i*pi * ((0:m-1)' + 1/3)/m );


% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(t,t,F');
contour(t,t,F',linspace(min(F(:)), max(F(:)),r), 'k');
colormap(parula(r-1));
caxis([min(F(:)) max(F(:))]);
axis image; axis off;

C = distinguishable_colors(length(Z));

tau = .01;
niter = 400;
for k=1:length(Z)
    z = Z(k);
    for i=1:niter
        z(end+1) = z(end) - tau*GradF(z(end));    
    end
    % plot curves
    plot(real(z), imag(z), 'LineWidth', 2, 'color', C(k,:));
end


saveas(gcf, [rep 'minmax-' num2str(round(10*a)) '.png'], 'png');

return;
z = z0;
for i=1:niter
    z(end+1) = z(end) - tau*GradF(z(end));    
end