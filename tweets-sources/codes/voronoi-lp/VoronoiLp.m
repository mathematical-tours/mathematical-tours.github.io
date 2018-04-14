%%
% Vornoi for Lp cost.

addpath('../toolbox/');
rep = MkResRep();

if not(exist('Y0'))
    % click and play
    Y0 = [];
    clf; hold on;
    while true
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 20);
        if button==3
            break;
        end
        Y0(:,end+1) = [a;b];
    end
end
k = size(Y0,2);


w = zeros(k,1);

q = 50;
plist = linspace(.3,5,q); 
for i=1:q
    t=(i-1)/(q-1);
    p = plist(i);
    DispVornoiLp(Y0,p);
    plot([0 1 1 0 0], [0 0 1 1 0], 'Color', [t, 0, 1-t], 'LineWidth', 5);
    drawnow;
    saveas(gcf, [rep 'vornoi-lp-' num2str(k) '-' znum2str(i,2) '.png'], 'png');
end

% AutoCrop(rep, ['vornoi-lp-' num2str(k)]); 

% convert iter-*.png lloyd.gif