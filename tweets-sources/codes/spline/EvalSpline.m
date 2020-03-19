function g = EvalSpline(f)

t = linspace(0,1,200);
switch length(f)
    case 2
        g = (1-t)*f(1) + t*f(2);
    case 3
        g = (1-t).^2*f(1) + 2*t.*(1-t)*f(2) + t.^2*f(3);
    case 4
        g = (1-t).^3*f(1) + 3*(1-t).^2.*t*f(2) + 3*(1-t).*t.^2*f(3)  + t.^3*f(4);
    otherwise
        error('Not implemented');
end

end