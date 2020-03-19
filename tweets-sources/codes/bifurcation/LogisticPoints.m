
addpath('../toolbox/');
rep = MkResRep();

n = 512; % #points
m = 128; % #param values
n1 = 300; % resolution for display


mymap = @(X).5 - abs(X-.5);
r = linspace(1.01,2, m); 

mymap = @(X)X .* (1-X);
r = linspace(2.4,4, m); 


x = linspace(0.01,.99,n);


[R,X] = meshgrid(r,x);

q = 100;
% precompute limit points 
for it=1:q    
   
        s = ones(n*m,1)*10; % size
    clf; scatter( R(:), X(:), s, 'filled' ); %  c(:,3),s, c, 'filled' );
    drawnow;
    
    tau = 0;
    X = tau*X + (1-tau) * R .* mymap(X);
    
    % drawnow;
    % imwrite(rescale(1-Xh), [rep 'anim-' znum2str(it,3) '.png']);
 
end

return;

clf; hold on;
for i=1:m
    s = (i-1)/(m-1);
    plot(r(i), X(i,:), '.', 'Color', [s 0 1-s]);
end