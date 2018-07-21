%%
% Mean Curvature motion evolution.

addpath('../toolbox/');
if not(exist('test'))
    test = 0;
end
test = test+1;
rep = MkResRep(num2str(test));

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


% click and play
if 1 % not(exist('gamma0'))
clf; hold on;
z0 = [];
while true
    axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    if button==3
        break;
    end
    z0(end+1) = a+1i*b;
end
gamma0 = z0(:);
n = length(z0);
end

% # points
p = 256*2;

% helpers
curvabs = @(gamma)[0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
resample1 = @(gamma,d)interp1(d/d(end),gamma,(0:p-1)'/p, 'linear');
resample = @(gamma)resample1( [gamma;gamma(1)], curvabs( [gamma;gamma(1)] ) );


gamma1 = resample(gamma0);
clf;
h = plot(gamma1([1:end 1]), 'k', 'LineWidth', 2);
axis('tight'); axis equal;  axis('off');

BwdDiff = @(c)c - c([end 1:end-1]);
FwdDiff = @(c)c([2:end 1]) - c;
dotp = @(c1,c2)real(c1.*conj(c2));
normalize = @(v)v./max(abs(v),eps);
tangent = @(gamma)normalize( FwdDiff(gamma) );
normal = @(gamma)-1i*tangent(gamma);
normalC = @(gamma)BwdDiff(tangent(gamma)) ./ abs( FwdDiff(gamma) );

dt = 0.001 / 500;
Tmax = 4 / 100;
niter = round(Tmax/dt);

gamma = gamma1;
q = 70; % frames
displist = round(linspace(1,niter,q));
k = 1;
for i=1:niter
    gamma = resample( gamma + dt * normalC(gamma) );
    if i==displist(k)
        t = (k-1)/(length(displist)-1);
        % display
        clf; hold on;
        h = plot(gamma([1:end 1]), 'r', 'color', [t 0 1-t], 'LineWidth', 2);
        if i==1 || i==niter
            set(h, 'LineWidth', 2);
        end
        axis('equal'); axis([-1 1 -1 1]); axis('off');
        drawnow;
        saveas(gcf, [rep 'mcm-' znum2str(k,2) '.png'], 'png');
        k = k+1;
    end
end


% AutoCrop(rep, 'mcm-');