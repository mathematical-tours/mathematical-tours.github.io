function [P,f,g,Err] = sinkhorn(C,a,b,epsilon,options)

a = a(:); b = b(:)';
n = length(a);

mina = @(H,epsilon)-epsilon*log( sum(a .* exp(-H/epsilon),1) );
minb = @(H,epsilon)-epsilon*log( sum(b .* exp(-H/epsilon),2) );
mina = @(H,epsilon)mina(H-min(H,[],1),epsilon) + min(H,[],1);
minb = @(H,epsilon)minb(H-min(H,[],2),epsilon) + min(H,[],2);

niter = getoptions(options, 'niter', 1000);
tol = getoptions(options, 'tol', 0);


f = zeros(n,1);
Err = [];
for it=1:niter
    g = mina(C-f,epsilon);
    f = minb(C-g,epsilon);
    % generate the coupling
    P = a .* exp((f+g-C)/epsilon) .* b;
    % check conservation of mass
    Err(it) = norm(sum(P,1)-b,1); 
    if Err(end)<tol
        break
    end
end

end