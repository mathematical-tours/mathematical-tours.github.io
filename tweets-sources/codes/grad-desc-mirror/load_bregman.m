function [psi,psii,phi,D] = load_bregman(bregmode)

switch bregmode
    case 'kl'
        psi  = @(x)log(x);
        psii = @(x)exp(x);
        phi = @(x)x.*(log(x)-1);
    case 'eucl'
        psi  = @(x)x;
        psii = @(x)x;
        phi = @(x)x.^2/2;
    case 'log'
        psi  = @(x)-1./x;
        psii = @(x)-1./x;
        phi = @(x)-log(x);
    case 'sqrt'
        psi  = @(x)-1./(2*sqrt(x));  
        psii = @(x)1./(2*x).^2;
        phi = @(x)-sqrt(x);
end
D = @(x,y)phi(x)-phi(y)-(x-y).*psi(y);

end