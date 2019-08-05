%%
% Interpolation using Shepard weights.

addpath('../toolbox/');
rep = MkResRep();

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

%% 
% In 1D.

% points and values
m = 6;  % #points
n = 256*4; % grid

rand('state', 123);
x = rand(m,1)+.5;
x = rescale(cumsum(x), .1, .9); 
f = rescale(rand(m,1));


X = linspace(0,1,n);
D = distmat(X(:)', x');


% same but movies

R = 40;
qlist = linspace(.2,5, R);
clf; hold on;
for i=1:R
 	clf; hold on;
    for j=1:i
        q = qlist(j);
        t=(j-1)/(R-1);
        %
        phi = @(r)1./(r+1e-10).^q;
        F = ( phi(D)*f ) ./ sum( phi(D), 2 );
        %
        X1 = [X(:);x(:)]; F1 = [F(:);f(:)]; [X1,I] = sort(X1); F1 = F1(I);
        hh = plot(X1,F1, 'LineWidth', 2, 'Color', [t 0 1-t]);
        if j<i
            hh.Color(4)=.2;
        end
    end
	plot(x,f, '.k', 'MarkerSize', 25);
    axis tight;
    box on;
    SetAR(2/3);
    axis([0,1,-.05,1.05]);
    set(gca, 'XTick', [], 'YTick', []); drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,3) '.png'], 'png');
end
% AutoCrop(repm, ['shepard-1d-']); 


