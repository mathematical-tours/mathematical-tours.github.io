%%
% display in colors of complex polynomials

addpath('../toolbox/');
rep = MkResRep();

% grid
r = 1.5;
n = 501;
t = linspace(-r,r,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;

name = 'five';
name = 'nine';
name = 'custom';
name = 'cauchy';
name = 'rational';
name = 'rational1';
name = 'ratioclick';

P = []; Pc = []; Proots = []; Ppoles = [];
switch name
    case 'five'
        P =@(z)z.^5-1;
        Pc = [1 0 0 0 0 -1];
    case 'nine'
        P =@(z)z.^9-1;
        Pc = [1 0 0 0 0 0 0 0 0 -1];
        
    case 'cauchy'
        Proots = [1i -1i]; Ppoles = [1 -1];
        F = (Z-1).*(Z+1)./(1+Z.^2);
        F = (Z.^2-1).*(1+Z.^2);
        
    case 'rational'
        F = Z./(1-Z.^5);
        Pc = [1 0 0 0 0 -1];
        Ppoles = 0;
        
    case 'rational1'
        F = (1-Z.^5) ./ Z;        
        F = (1-Z.^4) ./ (1-Z.^5);
        F = Z.*(1+Z.^5) ./ (1-Z.^5);
        Ppoles = roots([1 0 0 0 0 -1]);
        Proots = [roots([1 0 0 0 0 1]); 0];
        
        
    case 'ratioclick'
        
        F = ones(n);
        clf; hold on;
        while true
            axis equal; axis([-r r -r r]);
            if length(Proots)>0
               PolyDisp(F,t);
            end
            [a,b,button] = ginput(1);
            plot(a,b, 'b.', 'MarkerSize', 15);
            if button==3
                break;
            end
            Proots(end+1) = a+1i*b;
            F = F .* (Z-(a+1i*b));
        end
        while true
            axis equal; axis([-r r -r r]);
            PolyDisp(F,t);
            [a,b,button] = ginput(1);
            plot(a,b, 'r.', 'MarkerSize', 15);
            if button==3
                break;
            end
            Ppoles(end+1) = a+1i*b;
            F = F ./ (Z-(a+1i*b));
        end
        
        
        
        
    case 'custom'
        F = ones(n);
        clf; hold on;
        while true
            axis equal; axis([-r r -r r]);
            if length(Proots)>0
               PolyDisp(F,Proots,t);
            end
            [a,b,button] = ginput(1);
            plot(a,b, '.', 'MarkerSize', 15);
            if button==3
                break;
            end
            Proots(end+1) = a+1i*b;
            F = F .* (Z-(a+1i*b));
        end
end

if not(isempty(Pc))
    Proots = roots(Pc);
end
if not(isempty(P))
    F = P(Z);
end

clf; 
PolyDisp(F,t);
% roots plot/poles
plot(real(Proots), imag(Proots), 'b.', 'MarkerSize', 30);
plot(real(Ppoles), imag(Ppoles), 'r.', 'MarkerSize', 30);
saveas(gcf, [rep name '.png']);

