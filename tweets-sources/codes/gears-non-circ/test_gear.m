
%%
% Load a primal gear.

n = 1024*4;  % # sampling points

name = 'engrenage-16';
name = 'discont';
name = 'discont-2';
name = 'circle';
name = 'ellipse-centered';
name = 'hexagon';
name = 'random-strong';
name = 'petal';
name = 'petal-8';
name = 'cardioid';
name = 'square';
name = 'random';
name = 'triangle';
name = 'ellipse-focal';

center = [0 0];


rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end

tooth.transition =.2;
tooth.nbr = 30;
tooth.height = .03;
smoothing = .01;
x = load_gear(name, n, center, tooth, smoothing);

%%
% Compute the dual gear.

[y,L,phi] = compute_dual_gear(x);

%%
% Plot the gears


clf;
hold on;
plot_gear(x, [0 0], 'k',0, [],0);
plot_gear(y, [L*1.1 0], 'k',1, [],0);
axis tight; axis equal;
set(gcf, 'Color', [1 1 1]);
saveas(gcf, [rep 'gears-' name '.eps']);

return;

%% 
% display & save gif file

gifname = [rep name '.gif'];
plot_gear_animation(x,y,phi,L, 100, gifname);
