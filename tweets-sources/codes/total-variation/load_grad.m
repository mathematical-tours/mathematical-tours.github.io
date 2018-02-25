function [grad,div,lapl] = load_grad(N, bound)

% load_grad - load grad and div implicit operators
%
%   [grad,div,lapl] = load_grad(N);
%
% Copyright (c) 2016 Gabriel Peyre

if nargin<2
    bound = 'per';
    % bound = 'sym';
end
switch bound
    case 'per'
        s = [2:N 1];
        t = [N 1:N-1];
    case 'sym'
        s = [2:N N];
        t = [1 1:N-1];
end

grad = @(f)cat(3, ...
        (f-f(s,:)), ...
        (f(:,s)-f(s,s)), ...
        (f-f(:,s)),...
        (f(s,:)-f(s,s))  );
div = @(z)  ...
    -(z(:,:,1)-z(t,:,1)) + ...
    -(z(:,t,2)-z(t,t,2)) + ...
    -(z(:,:,3)-z(:,t,3)) + ...
    -(z(t,:,4)-z(t,t,4));
lapl = @(f)div(grad(f));

end
