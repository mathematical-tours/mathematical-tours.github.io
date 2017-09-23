% gradient descent step size
tau = 1.9/norm(Phi)^2;
% iterative thresholding
fspars = y; % initialization
for i=1:niter
	fspars = perform_thresholding(fspars + tau*Phi'*(y-Phi*fspars), lambda*tau, 'soft');
end