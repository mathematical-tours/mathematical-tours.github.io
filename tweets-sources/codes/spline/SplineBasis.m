function B = SplineBasis(x,t,n,p)

% t are the knots

p = length(t);
t = [t(:);t(:)+1]; % avoid boundary problems

% compute the basis
B = [];
for i=1:p
    B(:,i) = (x>=t(i)) & (x<t(i+1));
end
for d=1:n
    B1 = B;
    for i=1:p
        j = mod(i,p)+1;
        B1(:,i) = ...
            mod( (x-t(i)) , 1)/(t(i+d)-t(i)) .* B(:,i) ...
            + mod( (t(i+d+1)-x), 1)/(t(i+d+1)-t(i+1)) .* B(:,j);
    end
    B = B1;
end


end
