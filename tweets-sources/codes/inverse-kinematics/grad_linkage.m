function [dt,dx,L] = grad_linkage(r,t,y)

n = length(r);
x = forward_pass(r,t);
% loss function, computation of gradient using backprop
L = 1/2*abs(x(n+1)-y)^2; % loss 
dx(n+1) = x(n+1)-y; % dL/dx_{n+1}
for i=n:-1:1
    % x(i+1) = x(i) + r(i)*exp(1i*t(i));
    % ==> x(i+1) = x(i) + r(i)*[cos(t(i)), sin(t(i))];
    % dL/dx_i = dL/x_{i+1}*dx_{i+1}/dx_i 
    % dL/dt_i = dL/x1_{i+1}*dx1_{i+1}/dt_i  + dL/x2_{i+1}*dx1_{i+1}/dt_i 
    dx(i) = dx(i+1);
    dt(i) = -real(dx(i+1)) * r(i) * sin(t(i)) + ...
             imag(dx(i+1)) * r(i) * cos(t(i));
end

dx = dx(:); dt = dt(:);

end

