function y = eval_spline(m, x)

% eval_spline - evalutate a spline at an arbitrary location.
%
% y = eval_spline(m, x);
%
%   Copyright (c) 2004 Gabriel Peyré


if 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % direct evaluation
    y = x .* 0;
    fm = factorial(m);
    
    for i=0:m+1
        t = x + (m+1)/2 - i;
        t = t.^m .* (t>0) / fm;
        y = y + (-1)^i * nchoosek(m+1, i) * t;
    end
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % recursive evaluation
    if m<=0
        y = (x<1/2) .* (x>=-1/2);
    else
        mm = (m+1)/2;
        y = (mm+x)/m .* eval_spline(m-1,x+1/2) + (mm-x)/m .* eval_spline(m-1,x-1/2);
    end
    
end