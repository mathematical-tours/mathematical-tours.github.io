%% RELATED https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations#cite_note-23

nS = 250;
nI = 200;

s = linspace(1e-3, 2, nS);
i = linspace(0, 2, nI);
[I,S] = meshgrid(i,s);

g = 1;

M = I + S - g*log(S);
Mrange = [0 5];

% q = 50;
% v = linspace(0,1,

% display quantized colormap
r = 15; % #levellines
clf; hold on;
imagesc(s,i,M');
contour(s,i,M',linspace(Mrange(1),Mrange(2),r), 'k');
colormap(parula(r-1));
caxis(Mrange);
axis image; 
box on;
% axis off;
ylabel('S');
xlabel('I');
drawnow;
saveas(gcf, ['sir-example.png']);