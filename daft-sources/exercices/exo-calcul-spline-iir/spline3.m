function y = spline3(xx)
% calcul le spline de degre 3.
y = zeros(size(xx));

for i=1:length(xx)   
    x = xx(i);
    if( abs(x)>=0 && abs(x)<=1 )
        y(i) = 2/3-abs(x).^2+abs(x).^3/2;
    elseif( abs(x)<= 2 )
        y(i) = (2-abs(x))^3/6;
    else
        y(i) = 0;
    end
end