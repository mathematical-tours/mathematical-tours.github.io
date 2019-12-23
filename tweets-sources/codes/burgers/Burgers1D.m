%%
% Burgers equation in 1D
% df/dt = -d/dx( f^2 )

name = 'sine';
name = 'gauss';
name = 'box';

addpath('../toolbox/');
rep = MkResRep(name);

n = 256*2;

% initialization
t = (0:n-1)'/n;





switch name
    case 'gauss'
        % gaussian
        s = .06;
        f0 = exp(-(t-.2).^2/(2*s^2));
        Tmax = 2;
    case 'sine'
        f0 = sin(2*pi*t);
        Tmax = 1.5;
    case 'box';
        f0 = abs(t-.15)<.1;
        f0 = abs(t-.3)<.05;
        Tmax = 2.3;
        xi = .01/5; % diffusivity
end


xi = .1; % diffusivity
Tmax = .2;




% high diffusion
xi = 2; 
Tmax = .03;

% low diffusion
xi = .01/7; 
Tmax = 2.3;


% medium diffusion
xi = .015; 
Tmax = 2;



tau = 1/xi * .1/n^2;

niter = round(Tmax/tau);

dx = @(f)n * ( f([2:end 1])-f([end 1:end-1]) )/2;
ddx = @(f)n^2 * ( f([2:end 1])+f([end 1:end-1])-2*f ) / 2;

dx = @(f)n * ( f([2:end end])-f([1 1:end-1]) )/2;
ddx = @(f)n^2 * ( f([2:end end])+f([1 1:end-1])-2*f ) / 2;



q = 40;
ndisp = round(linspace(1,niter,q));
kdisp = 1;

f = f0;
% f = f-mean(f);
F = [];
for it=1:niter
    f = f + tau*( -dx( f.^2 )/2 + xi*ddx(f) );
    % f = f + tau*(  xi*ddx(f) );
    %
    if ndisp(kdisp)==it
        F(:,end+1) = f;
        clf; hold on;
        for jt=1:size(F,2)
            s = (jt-1)/(q-1);
            plot(t,F(:,jt), 'color', [s 0 1-s]*.3 + [1 1 1]*.7, 'LineWidth', 2);
        end
        plot(t,F(:,jt), 'color', [s 0 1-s], 'LineWidth', 3);
        axis tight;
        set(gca, 'Xtick', [], 'Ytick', [], 'PlotBoxAspectRatio', [1 1/2 1]); box on;
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(kdisp,2) '.png'] );
        kdisp = kdisp+1;
    end
    
end
% AutoCrop(rep, 'anim');