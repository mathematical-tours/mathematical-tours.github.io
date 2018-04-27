function [w,b,pi] = SolveSVM(X,lambda)

% generate class data
Z = [X{1},X{2}]';
n = size(Z,1);
y = ones(n,1); y(size(X{1},2)+1:end) = -1;

% max 1/n sum max(0,1-y*(x_i*w-b) + lambda*|w|^2
% 1/n * sum( max(0, 1-y.*(Z*w-b)) ) + lambda*norm(w)^2

svm_mode = 'hard';
svm_mode = 'soft';

if lambda==0
    %% HARD MARGIN SVM %%
    cvx_begin quiet
    variables w(2,1) b;
    dual variables pi;
    minimize( .5*w'*w );
    subject to
    pi : y .* ( Z*w + b) >= 1; %
    cvx_end
else
    %% SOFT MARGIN %%%
    cvx_begin quiet
    variables w(2,1) zeta(n,1) b;
    dual variables pi;
    minimize( 1/n * sum(zeta) + lambda * .5*w'*w );
    subject to
    pi : y .* ( Z*w + b) >= 1-zeta; %
    zeta >= 0;
    cvx_end
end

end

