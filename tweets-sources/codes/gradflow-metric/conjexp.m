function q = conjexp(p)

% conjexp
%
%   1/p+1/q=1


q = p/(p-1);
if p==1
    q = Inf;
end
if p==Inf
    q = p;
end

end