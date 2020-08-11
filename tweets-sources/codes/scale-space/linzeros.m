function [ti,gi] = linzeros(t,f, g)

t = t(:); f = f(:); g = g(:);

I = find( sign(f([2:end])) .* sign(f([1:end-1]))  < 0 );

t1 = t(I); t2 = t(I+1);
f1 = f(I); f2 = f(I+1);
g1 = g(I); g2 = g(I+1);

% f1 + (f2-f1)/(t2-t1) * (t-t1) = 0;
% t = t1 + -f1 * (t2-t1)/(f2-f1)

ti = t1 + -f1 .* (t2-t1)./(f2-f1);
gi = g1 + (g2-g1)./(t2-t1) .* (ti-t1);


end