function v = randdisc(p,Q)

% randdisc - discrete random number generator.
%
%   v = randdisc(p,Q)
%
%   Copyright (c) 2015 Gabriel Peyre

normalize = @(p)p./repmat(sum(p,1),[size(p,1) 1]);

N = size(p,1); % size of distributions
S = size(p,2); % number of distributions
p = normalize(p+max(p(:))*1e-9);

a = [zeros(1,S); cumsum(p,1)];
b = (0:N)' * ones(1,S);
u = rand(Q,S);
for i=1:S
    v(:,i) = interp1(a(:,i),b(:,i),u(:,i),'next');
end

end