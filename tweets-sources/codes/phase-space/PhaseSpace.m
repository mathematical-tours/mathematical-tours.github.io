%%
% Display evolution in phase space.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

name = 'pendulum';
name = 'oscilator';

q = 1; % number of parameter to tests
m = 8; % number of initial conditions
r = 50; % number of frames of the anim

% for savings
for k=1:m
    [~,~] = mkdir([rep '/' num2str(k) '/']);
end

% 3 / 6 / 10 / 13

switch name
    case 'oscilator'
        v = linspace(-1,1,20);
        x = linspace(-1,1,20);
        % param list
        %  0 : pure circles
        %  smaller than 1 : oscillation
        zeta = .8; 
        % initial points
        x0_list = linspace(.2,.98,m)'; 
        v0_list = zeros(m,1); 
        %
        x0_list = linspace(.1,.7,m)'; 
        v0_list = linspace(.1,.7,m)'; 
        %
        omega = 1;
        tmax = 10;
        
    case 'pendulum'
        v = linspace(-2,2,20);
        x = linspace(-4,4,24);
        % param list
        zeta = + .5;
        % initial points
        x0_list = linspace(0,-4,m+1);  % 0 * ones(m,1);
        v0_list = linspace(0,1.6,m+1);
        x0_list(1) = []; v0_list(1) = [];
        
        
        x0_list = linspace(4,-4,m);  % 0 * ones(m,1);
        v0_list = linspace(-1.6,1.6,m);
        tmax = 100;
        
        
end

[V,X] = meshgrid(v,x);
nt = 6000; % time steps
t = linspace(0,tmax,nt)';



Col = distinguishable_colors(m);

% vector field
switch name
    case 'oscilator'
        Ax = V;
        Av = -2*zeta*omega*V - omega^2*X;
    case 'pendulum'
        Ax = V;
        Av = -zeta*sin(X);
end
Ax1 = Ax ./ sqrt(Ax.^2+Av.^2); % normalize
Av1 = Av ./ sqrt(Ax.^2+Av.^2);
% integrate in time
xt = zeros(nt,m);
vt = zeros(nt,m);
for i=1:m
    x0 = x0_list(i);
    v0 = v0_list(i);
    switch name
        case 'oscilator'
            %
            % some trajectories
            m1 = -zeta*omega + omega*sqrt(zeta^2-1);
            m2 = -zeta*omega - omega*sqrt(zeta^2-1);
            %
            a = inv([1 1; m1 m2])*[x0; v0];
            %
            xt(:,i) = real( a(1)*exp(t*m1) + a(2)*exp(t*m2) );
            vt(:,i) = real( a(1)*m1*exp(t*m1) + a(2)*m2*exp(t*m2) );
        case 'pendulum'
            % time integration
            tau = .002;
            xt(1,i) = x0; vt(1,i) = v0;
            for k=1:nt-1
                [xt(k+1,i),vt(k+1,i)] = deal( ...
                    xt(k,i)+tau*vt(k,i), vt(k,i) + tau*( -zeta*sin(xt(k,i)) ) );
            end
    end
    
end

% display evolution
idisp = round(linspace(1,nt,r));
for it=1:r
    imax = idisp(it);
    % vector field
    clf; hold on;
    quiver(x, v, Ax1', Av1', ...
        'k', 'filled', 'LineWidth', 1, 'AutoScaleFactor', .7);
    %
    for j=1:m
        col = [(j-1)/(m-1), 0, 1-(j-1)/(m-1)];
        plot( xt(1:imax,j), vt(1:imax,j), 'color',col, 'LineWidth', 3 );
        plot( xt(imax,j), vt(imax,j), '.', 'color',col, 'MarkerSize', 20 );
    end
    axis equal; axis off;
    axis([min(x) max(x) min(v) max(v)]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);
    % pendulum
    % display pendulum
    for j=[2 5 8]% 1:m %
        col = [(j-1)/(m-1), 0, 1-(j-1)/(m-1)];
        switch name
            case 'pendulum'
                % [~,~] = mkdir([rep '/' num2str(k) '/']);
                clf; hold on;
                plot([0 exp(1i * (pi/2-xt(imax,j)))], '.-', 'color',col, 'LineWidth', 3, 'MarkerSize', 20);
                plot([exp(1i * (pi/2-xt(imax,j)))], '.', 'color',col, 'LineWidth', 3, 'MarkerSize', 50);
                axis ij; axis equal; 
                axis([-1 1 -1 1]*1.1);
            case 'oscilator'
                u = linspace(0,5,2000);
                clf; hold on;
                plot(xt(imax,j)*u/max(u),(abs(2*mod(u-.25,1)-1)-.5)*.5, '.-', 'color',col, 'LineWidth', 3, 'MarkerSize', 20);
                plot(xt(imax,j),0, '.', 'color',col, 'LineWidth', 3, 'MarkerSize', 50);
                axis([-1 1 -.3 .3]*1.1);
                set(gca, 'PlotBoxAspectRatio', [1 .2 1])
        end
        box on;  set(gca, 'XTick', [], 'YTick', []);
        saveas(gcf, [rep '/' num2str(j) '/' 'anim-' znum2str(it,3) '.png']);
        % display curves
        clf; hold on;
        plot(t, xt(:,j), 'k-');
        plot(t(1:imax), xt(1:imax,j), '-', 'color', col, 'LineWidth', 3);
        plot(t(imax), xt(imax,j), '.', 'color', col, 'MarkerSize', 25);
        axis([0 max(t) -1 1]);
        box on;
        set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/2 1]);
        saveas(gcf, [rep '/' num2str(j) '/' 'curve-' znum2str(it,3) '.png']);
        
    end
end

