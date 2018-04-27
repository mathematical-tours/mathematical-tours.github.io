%%
% Display of SVM classification method.

addpath('../toolbox/');
rep = MkResRep();

% SVM with more and more points

lambda = 0.001;
name = 'soft1';

ms = 10;
col = {[0 0 1] [1 0 0]};
r = .3;
LL = @(c)r*c + (1-r)*[1 1 1]; % lighten a color


X = {[] []}; Xsvg = {};
iter = 0; 
clf; hold on; axis equal;
axis([0 1 0 1]); axis on; box on; set(gca, 'XTick', [], 'YTick', []);     
while true
	[a,b,button] = ginput(1);
    i = 1;
    if button==3
        i = 2;
    end
    if button==2
        break;
    end
    X{i}(:,end+1) = [a;b];
    %% display
    clf; hold on; axis equal;
    axis([0 1 0 1]); axis on; box on; set(gca, 'XTick', [], 'YTick', []); 
    if size(X{1},2)>0 && size(X{2},2)>0
        [w,b,pi] = SolveSVM(X,lambda);
        PlotSVM(X,w,b,pi);
    else
        for i=1:2
            if not(isempty(X{i}))
               	plot(X{i}(1,:), X{i}(2,:), 'o', 'MarkerFaceColor', LL(col{i}), 'MarkerSize', ms, 'MarkerEdgeColor', col{i});
            end
        end
    end
    %% save
	Xsvg{end+1} = {X};
    iter = iter+1;
    saveas(gcf, [rep name '-' znum2str(iter,2) '.png']);
end



lambda_list = [.01 0];
names = {'soft2' 'hard'};
for it=1:length(lambda_list)
    lambda = lambda_list(it);
    for iter=1:length(Xsvg)
        X = Xsvg{iter}{1};
        
        clf; hold on; axis equal;
        axis([0 1 0 1]); axis on; box on; set(gca, 'XTick', [], 'YTick', []);
        if size(X{1},2)>0 && size(X{2},2)>0
            [w,b,pi] = SolveSVM(X,lambda);
            PlotSVM(X,w,b,pi);
        else
            for i=1:2
                if not(isempty(X{i}))
                    plot(X{i}(1,:), X{i}(2,:), 'o', 'MarkerFaceColor', LL(col{i}), 'MarkerSize', ms, 'MarkerEdgeColor', col{i});
                end
            end
        end
        saveas(gcf, [rep names{it} '-' znum2str(iter,2) '.png']);
        
    end
end



names = {'soft1' 'soft2' 'hard'};
for it=1:length(names)
    AutoCrop(rep, [names{it} '-']);
end


% convert -delay 1x3 ....

