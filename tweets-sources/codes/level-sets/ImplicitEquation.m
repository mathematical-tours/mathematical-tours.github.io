function F = ImplicitEquation(x,y,z,name)


switch name
    case 'ball'
        F = sqrt( x.^2 + y.^2 + z.^2 ) - .8;
    case 'torus'
        R = .7; r = .25;
        F = (x.^2+y.^2+z.^2 + R^2-r^2).^2 - 4*R^2*(x.^2+y.^2);
        F = abs(x.^2+y.^2+z.^2 + R^2-r^2) - 2*R*sqrt(x.^2+y.^2);
    case 'torus2'
        x = x*1.2+1;
        z = z/2;
        F = (x.*(x-1).^2.*(x-2) + y.^2).^2 + z.^2 - .1.^2;
    case 'trumpet'
        a = .5;
        F = (x.^2 + y.^2).*z.^2 - a^4;
    case 'pball'
        p = .7;
        F = ( abs(x).^p + abs(y).^p + abs(z).^p ) - 1;
end