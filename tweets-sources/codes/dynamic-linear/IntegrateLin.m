function y = IntegrateLin(A,x0,t)

n = length(t); 
m = length(x0);
x = [real(x0),imag(x0)];

if 0
    %% SLOW %%
    y = zeros(2,m,n);
    for i=1:n
        y(:,:,i) = expm(t(i)*A)*x';
    end
else
    %% FAST %%
    % A = U*S*inv(U)
    % exp(t*A) = U*diag(exp(t*S))*inv(U)
    [U,S] = eig(A);
    S = diag(S);
    M = exp( S(:) .* reshape(t(:), [1 1 n] ) ) .* ( inv(U)*x' );
    y = reshape( U*reshape(M, [2 m*n]), [2 m n] );
end
y = transpose( squeeze( y(1,:,:) + 1i*y(2,:,:) ) );



end
