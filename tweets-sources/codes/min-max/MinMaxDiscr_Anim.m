%%
% Test for min-max games optimization.

addpath('../toolbox/');
rep = '../results/min-max/';
[~,~] = mkdir(rep);

q = 501;
t = linspace(-1,1,q);
[Y,X] = meshgrid(t,t);

desc_mode = 'implicit';
desc_mode = 'extra';
desc_mode = 'explicit';

averaging = 0;


F = X.*Y + .1*X.^2/2;
GradF = @(z)(imag(z)+.1*real(z))  - 1i * real(z);

F = X.*Y;
GradF = @(z)imag(z)  - 1i * real(z);



z0 = .5*exp( 2i*pi/3 );

%% Explicit Euler %%

% step sizes
tau_list = [.001 .1 .2 .3];
tau_list = [.1 .2];
C = distinguishable_colors(length(tau_list));


tau = .1;
tau = .2;

niter = 50;
q = 50; % for anim
ndisp = round(linspace(1,niter,q));
kdisp = 1;

col = [1 0 0];

C = distinguishable_colors(5);

C = [1 0 0; 0 1 0; 0  0 1];


       
       
    
% display quantized colormap
z = [z0 z0 z0]; za = z;
for it=1:niter
    % 'explicit'
    z(it+1,1) = z(it,1) - tau*GradF(z(it,1));
    % 'extra'
    z1 = z(it,2) - tau*GradF(z(it,2));
    z(it+1,2) = z(it,2) - tau*GradF(z1);
    % . 'implicit'
    u = [1, tau; -tau, 1]\[real(z(it,3)); imag(z(it,3))];
    z(it+1,3) = u(1) + 1i*u(2);
    if averaging==0
        za(it+1,:) = z(it+1,:);
    else
        za(it+1,:) = mean(z);
    end
    
    if ndisp(kdisp)==it
        % display
        r = 15; % #levellines
        clf; hold on;
        imagesc(t,t,F');
        contour(t,t,F',linspace(min(F(:)), max(F(:)),r), 'k');
        colormap(parula(r-1));
        caxis([min(F(:)) max(F(:))]);
        [x,y] = meshgrid(-1:.15:1,-1:.15:1);
        quiver(x,y,-y,x,'k');
        axis image; axis off;
        % plot curves
        for k=[1 3]
            col = C(k,:);
            plot(real(za(:,k)), imag(za(:,k)), '.-', 'LineWidth', 2, 'color', col, 'MarkerSize', 15);
            plot(real(za(end,k)), imag(za(end,k)), '.', 'LineWidth', 2, 'color', col, 'MarkerSize', 30);
        end
        axis equal; axis([-1 1 -1 1]);
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
        kdisp = kdisp+1;
    end
    
end


% saveas(gcf, [rep 'discr-' desc_mode '.png'], 'png');
