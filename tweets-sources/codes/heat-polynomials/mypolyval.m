function F = mypolyval(XY,r)

F = ones(size(XY));
for i=1:length(r)
    F = F .* (XY-r(i));
end

end