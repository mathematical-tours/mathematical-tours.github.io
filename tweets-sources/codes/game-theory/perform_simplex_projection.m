function y = perform_simplex_projection(x,rho)
% perform_simplex_projection - compute the projection on the simplex of size rho
%
%   y = perform_simplex_projection(x,rho);
%
%   x is the projection of y on the set {a \ a >= 0 sum_i a_i = rho },
%
%   Copyright (c) 2013 Jalal Fadili
    

[n,d] = size(x);

if rho<=0
    error('rho should be > 0');
end

if ~isvector(rho) | numel(rho)~=n
    error('rho must be a vector of size n, x is n x d');
end

rho = repmat(rho(:),[1 d]);
	
xs = sort(x,2,'descend');
xtmp = (cumsum(xs,2)-rho)*diag(sparse(1./(1:d)));
y = max(bsxfun(@minus,x,xtmp(sub2ind([n,d],(1:n)',sum(xs>xtmp,2)))),0);

end


