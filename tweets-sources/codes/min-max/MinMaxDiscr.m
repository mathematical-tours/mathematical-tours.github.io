%%
% Test for min-max games optimization.

addpath('../toolbox/');
rep = '../results/min-max/';
[~,~] = mkdir(rep);

q = 501;
t = linspace(-1,1,q);
[Y,X] = meshgrid(t,t);

desc_mode = 'implicit';
desc_mode = 'explicit';
desc_mode = 'extra';

averaging = 0;


F = X.*Y;
GradF = @(z)imag(z)  - 1i * real(z);



z0 = .4*exp( 2i*pi/3 );

%% Explicit Euler %%

% step sizes
tau_list = [.001 .1 .2 .3];
tau_list = [.1 .2];
C = distinguishable_colors(length(tau_list));

% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(t,t,F');
contour(t,t,F',linspace(min(F(:)), max(F(:)),r), 'k');
colormap(parula(r-1));
caxis([min(F(:)) max(F(:))]);
axis image; axis off;
for i=1:length(tau_list)
    tau = tau_list(i);
    niter = round(3*pi/tau);
    niter = 200;
    z = z0; za = z0;
    for it=1:niter
        switch desc_mode
            case 'explicit'
                z(end+1) = z(end) - tau*GradF(z(end));
            case 'extra'
                z1 = z(end) - tau*GradF(z(end));
                z(end+1) = z(end) - tau*GradF(z1);
            case 'implicit'
                u = [1, tau; -tau, 1]\[real(z(end)); imag(z(end))];
                z(end+1) = u(1) + 1i*u(2);
        end
        if averaging==0;
            za(end+1) = z(end);
        else
            za(end+1) = mean(z);
        end
    end
    % plot curves
    plot(real(za), imag(za), '.-', 'LineWidth', 2, 'color', C(i,:), 'MarkerSize', 20);
end
axis([-1 1 -1 1])
% saveas(gcf, [rep 'discr-' desc_mode '.png'], 'png');
