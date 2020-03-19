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
while true
    
    clf; hold on;    
    if not(isempty(X))
        % optimize weights, and predict
        % 'rbf' 'polynomial'
        w  = fitcsvm(X,Y,'KernelFunction','polynomial', ...
            'KernelScale', .1, ...
            'PolynomialOrder', 3, ...
            'Prior', 'uniform');
        Y1 = reshape(predict(w,X1), p);
        clf; hold on;
        imagesc(s,t,Y1');
        scatter(X(w.IsSupportVector,1),X(w.IsSupportVector,2),'ok');
        scatter(X(Y==1,1),X(Y==1,2),'.b');
        scatter(X(Y==-1,1),X(Y==-1,2),'.r');        
    end
    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    %
    [a,b,button] = ginput(1);
    if button==3
        Y(end+1) = -1;
    else
        Y(end+1) = 1;
    end
    X(end+1,:) = [a b];
    
end
