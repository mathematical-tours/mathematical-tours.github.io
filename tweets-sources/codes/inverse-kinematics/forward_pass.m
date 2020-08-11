function x = forward_pass(r,t)

n = length(r);
% forward anim
x = 0;
for i=1:n
    x(i+1) = x(i) + r(i)*exp(1i*t(i));
end

end