% Initialization of the dual variables.
G = zeros(n,n,2);
for i=1:niter
    % Gradient of the energy
    dG = grad( div(G) - f/lambda );
    % gradient descent
    G = G + tau*dG;
    % Projection on Linfty constraints
    d = repmat( sqrt(sum(G.^2,3)), [1 1 2] ); % norm of the vectors
    G = G ./ max(d,ones(n,n,2));
end
% Final solution from the dual variables.
ftv = f-lambda*div(G);