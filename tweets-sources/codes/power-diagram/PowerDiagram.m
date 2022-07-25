%%
% test for k-means on a dense lattice.

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


w = .1 * (1:k)/k;
w = w(randperm(length(w)));

q = 50; 
for i=1:q
    t=(i-1)/(q-1);
    DisplayDiagram(Y0,t*w,(1-t)*1/2 + 1.5*t*w/max(w));
    drawnow;
    % saveas(gcf, [rep 'power-diag-' num2str(k) '-' znum2str(i,2) '.png'], 'png');
end

% AutoCrop(rep, ['power-diag-' num2str(k)]); 


% AutoCrop(repsvg, 'iter-')
% convert iter-*.png lloyd.gif