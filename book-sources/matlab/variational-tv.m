% step size
tau = epsilon/5;
% number of iterations to reach T.
niter = round(T/tau);
% initial solution for time t=0.
ftv = f;
for i=1:niter
	% gradient and norm of gradient
    Gr = grad(ftv); d = sqrt(sum3(Gr.^2,3));
    G = -div( Gr ./ repmat( sqrt( epsilon^2 + d.^2 ) , [1 1 2]) );
    % descent step
	ftv = ftv - tau*G;
end