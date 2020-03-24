%%
% Test for kernelize svm.


addpath('../toolbox/');
addpath('./svm/');
rep = MkResRep();

% grid for prediction
p = [200 200];
s = linspace(0,1,p(1));
t = linspace(0,1,p(2));
[T,S] = meshgrid(t,s);
X1 = [S(:),T(:)];


X = []; Y = [];
it = 0;
while true
    
    clf; hold on;    
    if not(isempty(X))
        % optimize weights, and predict
        % 'rbf' 'polynomial'
        w  = fitcsvm(X,Y,'KernelFunction','rbf', ...
            'KernelScale', .2, ...
            'Prior', 'uniform');
            % 'PolynomialOrder', 3, ...
        Y1 = reshape(predict(w,X1), p);
        clf; hold on;
        imagesc(s,t,Y1');
        scatter(X(w.IsSupportVector,1),X(w.IsSupportVector,2),'ok');
        scatter(X(Y==1,1),X(Y==1,2),'.b');
        scatter(X(Y==-1,1),X(Y==-1,2),'.r');        
    end
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    if it>=1
        saveas(gcf, [rep  'anim-' znum2str(it,3) '.png'] );
    end
    %
    [a,b,button] = ginput(1);
    if button==3
        Y(end+1) = -1;
    else
        Y(end+1) = 1;
    end
    X(end+1,:) = [a b];
    it = it+1;
    
end
