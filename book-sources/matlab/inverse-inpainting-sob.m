% step size 
tau = 1/5;
% initialization
fsob = y;
for i=1:niter
	% laplacien
    L = div( grad(fsob) );
    % gradient descent
    fsob = fsob + tau * L;
    % project on constraints
    fsob(mask==0) = y(mask==0);
end
