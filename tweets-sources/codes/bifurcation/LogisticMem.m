
addpath('../toolbox/');
rep = MkResRep();

n = 512; % #points
m = 128*4; % #param values
n1 = 300; % resolution for display


mymap = @(X).5 - abs(X-.5);
r = linspace(1.01,2, m); 

mymap = @(X)X .* (1-X);
r = linspace(2.4,4, m); 


x = linspace(0.01,.99,n);


[R,X] = meshgrid(r,x);

X = ones(1,m)*.5;

niter = 1000;
q = 60;
ndisp = unique( 1 + round((niter-1)*linspace(0,1,q).^3) );
q = length(ndisp);
mem = 300;

kdisp = 1;
% precompute limit points 
for it=1:niter  
   
	% compute histogram of appearance
    if it==1
        Xh = hist([X;X],linspace(0,1,n1));
    else
        Xh = hist(X,linspace(0,1,n1));        
    end
    Xh = Xh ./ max(Xh);
    
    if it==ndisp(kdisp)
    clf;
    imageplot(1-Xh);
    % imageplot(X);
    axis tight;
    drawnow;
    
	% drawnow;
    imwrite(rescale(1-Xh), [rep 'anim-' znum2str(kdisp,3) '.png']);
    kdisp = kdisp+1;
    end
    
    tau = 0;
    X(end+1,:) = tau*X(end,:) + (1-tau) * r .* mymap(X(end,:));
    X = X(max(size(X,1)-mem,1):end,:);
end

return;

clf; hold on;
for i=1:m
    s = (i-1)/(m-1);
    plot(r(i), X(i,:), '.', 'Color', [s 0 1-s]);
end