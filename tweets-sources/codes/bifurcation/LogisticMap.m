
addpath('../toolbox/');
rep = MkResRep();

n = 512*10; % #points
n1 = 300; % resolution for display
m = 128*4; % #param values

mymap = @(X)X .* (1-X);
r = linspace(2.4,4, m); 

mymap = @(X).5 - abs(X-.5);
r = linspace(1.01,2, m); 


x = linspace(0,1,n);


[R,X] = meshgrid(r,x);

q = 40;

for it=1:q
    
    tau = 0;
    X = tau*X + (1-tau) * R .* mymap(X);
   
    % compute histogram of appearance
    Xh = hist(X,linspace(0,1,n1));
    Xh = Xh ./ max(Xh);
    
    clf;
    imageplot(1-Xh);
    % imageplot(X);
    axis tight;
    drawnow;
    
	clf;
    imageplot(X);
    % imageplot(X);
    axis square;
    % drawnow;
    
    % imwrite(rescale(1-Xh), [rep 'anim-' znum2str(it,3) '.png']);
    
    
end