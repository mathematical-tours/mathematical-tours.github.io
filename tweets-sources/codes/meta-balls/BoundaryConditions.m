function [c,s] = BoundaryConditions(c,s, BB)

% reflexion boundary condition on a box

BBa = BB(1:2:end);
BBb = BB(2:2:end);

for k=1:3
    I = find(c(:,k)<BBa(k));
    c(I,k) = BBa(k);
    s(I,k) = -s(I,k);
    %
    I = find(c(:,k)>BBb(k));
    c(I,k) = BBb(k);
    s(I,k) = -s(I,k);
end

end