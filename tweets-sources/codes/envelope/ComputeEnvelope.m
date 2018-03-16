function h = ComputeEveloppe(c,cn, bound)

if nargin<3
    bound = 'per';
end

n = length(c);

if strcmp(bound, 'per')
    nmax = n;
else
    nmax = n-1;
end

for i=1:nmax
    j = mod(i,n)+1;
    m = [real(cn(i)),-real(cn(j)); imag(cn(i)),-imag(cn(j))] \ ...
        [real(c(j)-c(i)); imag(c(j)-c(i))];
    h(i) = c(i)+m(1)*cn(i);
    %h1(i) = c(j)+m(2)*cn(j);
end
h = h(:);

end