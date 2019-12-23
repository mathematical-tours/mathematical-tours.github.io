%%
% Burgers equation in 1D
% df/dt = -d/dx( f^2 )

name = 'heat';
name = 'wave';

addpath('../toolbox/');
rep = MkResRep(name);

n = 256;

% initialization
t = linspace(0,1,n); % (0:n-1)'/n;






f0 = abs(t-.3)<.05;
for i=1:30
    f0 = (f0 + .5*f0([1 1:end-1]) + .5*f0([2:end end]))/2;
end


m = .3; s = .03;
f0 = exp( -(t-m).^2/(2*s^2) );

switch name
    case 'heat'
        % for heat
        Tmax = .07;
        tau = .9/n^2;
        niter = round(Tmax/tau);
    case 'wave'
        % for wave
        Tmax = .01;
        tau = .05/n^2;
        niter = 1800;
end


dx = @(f)n * ( f([2:end end])-f([1 1:end-1]) )/2;
ddx = @(f)n^2 * ( f([2:end end])+f([1 1:end-1])-2*f ) / 2;



q = 50;
ndisp = round(linspace(1,niter,q));
kdisp = 1;

f = f0;
f1 = f0;
F = [];
for it=1:niter
    switch name 
        case 'heat'
            f = f + tau* ddx(f);
        case 'wave'
            [f,f1] = deal( 2*f-f1+tau*ddx(f), f );
    end
    %
    if ndisp(kdisp)==it
        F(:,end+1) = f;
        clf; hold on;
        dr = 1;
        if 1
        for jt=1:size(F,2)
            s = (jt-1)/(q-1);
            
            % 1 <-> jt=size(F,2)
            % 0 <-> jt=size(F,2)-7            
            dr = 1 - (size(F,2)-jt-12)/7;
            dr = min(max(dr,0),1);
            dr = .3; 
            plot(t,F(:,jt), 'color', [s 0 1-s]*dr + [1 1 1]*(1-dr), 'LineWidth', 2);
            dr = dr*.7;
        end
        end
        s = (size(F,2)-1)/(q-1);
        plot(t,F(:,end), 'color', [s 0 1-s], 'LineWidth', 3);
        axis tight;
        axis([0 1 -.02 1]);
        set(gca, 'Xtick', [], 'Ytick', [], 'PlotBoxAspectRatio', [1 1/2 1]); box on;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(kdisp,2) '.png'] );
        kdisp = kdisp+1;
    end
    
end
% AutoCrop(rep, 'anim');