
%%
% Load a primal gear.

n = 1024*4;  % # sampling points

name = 'discont';
name = 'circle';
name = 'ellipse-centered';
name = 'hexagon';
name = 'petal';
name = 'cardioid';
name = 'square';
%
%
name = 'engrenage-16';
name = 'petal-8';
name = 'discont-2';
name = 'triangle';
name = 'ellipse-focal';
name = 'random';
name = 'random-strong';


addpath('../toolbox/');
rep = MkResRep(name);

center = [0 0];


tooth.transition =.2;
tooth.nbr = 30;
tooth.nbr = 0; % no tooth
tooth.height = .03;
smoothing = .01;
x = load_gear(name, n, center, tooth, smoothing);

%%
% Compute the dual gear.

[y,L,phi] = compute_dual_gear(x);

%% 
% display

q = 50; % number of time step

b = max(max(x),max(y));
a1 = max(x); a2 = L+max(y);

t = 2*pi* (0:n-1)'/n;
tx = linspace(0,2*pi,q+1)';
ty = linspace(0,2*pi,n+1)';
ty = interp1( ty, phi, tx);

for i=1:q
    clf;
    hold on;
    plot_gear(x, [0 0], 'b',0, t + tx(i));
    plot_gear(y, [L 0], 'r',1, t + ty(i));
    axis([-a1-.1 a2+1 -b-.1 b+.1]); axis on;
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'XTick', [], 'YTick', []); box on;
    drawnow;
    % saveas(gcf, [rep 'anim-' znum2str(i,2) '.png' ]);
end

% AutoCrop(rep,'anim');
    

return;

%% 
% display & save gif file

gifname = [rep name '.gif'];
plot_gear_animation(x,y,phi,L, 100, gifname);
