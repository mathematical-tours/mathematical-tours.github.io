%%
% Display evolution in phase space.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

name = 'oscilator';
name = 'pendulum';

q = 50;
m = 16; % number of initial conditions

switch name
    case 'oscilator'
        v = linspace(-1,1,20);
        x = linspace(-1,1,20);
        % param list
        zeta_list = linspace(0, 3, q);
        % initial points
        x0_list = .8 * cos(2*pi*(1:m)/m);
        v0_list = .8 * sin(2*pi*(1:m)/m);
        omega = 1;

    case 'pendulum'
        v = linspace(-2,2,20);
        x = linspace(-4,4,24);
        % param list
        zeta_list = linspace(0, 3, q);
        % initial points
        x0_list = 0 * ones(m,1);
        v0_list = linspace(-1.9,1.9,m);
        zeta_list = linspace(.01, 1, q);
end

[V,X] = meshgrid(v,x);

for it=1:q
    zeta = zeta_list(it);
    switch name
        case 'oscilator'
            Ax = V;
            Av = -2*zeta*omega*V - omega^2*X;
        case 'pendulum'
            Ax = V;
            Av = -zeta*sin(X);
    end
    %
    Ax1 = Ax ./ sqrt(Ax.^2+Av.^2);
    Av1 = Av ./ sqrt(Ax.^2+Av.^2);
    clf; hold on;
        quiver(x, v, Ax1', Av1', ...
        'k', 'filled', 'LineWidth', 1, 'AutoScaleFactor', .7);
    for i=1:m
        x0 = x0_list(i);
        v0 = v0_list(i);
        switch name
            case 'oscilator'
                %
                % some trajectories
                m1 = -zeta*omega + omega*sqrt(zeta^2-1);
                m2 = -zeta*omega - omega*sqrt(zeta^2-1);
                t = linspace(0,100,5000);
                %
                a = inv([1 1; m1 m2])*[x0; v0];
                %
                xt = a(1)*exp(t*m1) + a(2)*exp(t*m2);
                vt = a(1)*m1*exp(t*m1) + a(2)*m2*exp(t*m2);
                xt = real(xt); vt = real(vt);
            case 'pendulum'
                % time integration
                tau = .002; nt = 10000;
                xt = x0; vt = v0;
                for k=1:nt
                    [xt(end+1),vt(end+1)] = deal( ...
                        xt(end)+tau*vt(end), vt(end) + tau*( -zeta*sin(xt(end)) ) );
                end
        end

        %
        s = sin(pi*i/m)^2;
        s = (it-1)/(q-1);
        col = [s, 0, 1-s];
        %
        plot( xt, vt, 'color',col, 'LineWidth', 3 );
        plot( x0, v0, '.', 'color',col, 'MarkerSize', 20 );
        axis equal; axis off;
        axis([min(x) max(x) min(v) max(v)]);
    end
    drawnow;
    %clf;
    %plot(xt);
    %hold on;
    %drawnow;
end
