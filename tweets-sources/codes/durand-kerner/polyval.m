function p = polyval(x,z)

% p=prod_{i}(x-z_i)

p = 1;
for i=1:length(z)
    p = p*(x-z(i));
end

end