function T = load_gabor_kernel(n, KerType, options)

% load_gabor_kernel - load a 2-D spacial kernel for synthesis
%
% T = load_gabor_kernel(n, KerType, options);
%
%   Copyright (c) 2012 Gabriel Peyre


Lx = [0:n/2-1 -n/2:-1];
[Y,X] = meshgrid(Lx,Lx);


%%
% Speed enveloppe.

direction = @(rho,theta)rho*[cos(theta);sin(theta)];
gauss = @(u,s)exp(-1/2*( u ./ max(abs(s),1e-10) ).^2 );


%%
% spacial enveloppe.

equalize = @(a)a/norm(a(:));
switch KerType
    case 'gabor'
        rho = getoptions(options, 'rho', 0, 1);
        theta = getoptions(options, 'theta', 0, 1);
        sigmax = getoptions(options, 'sigmax', 0, 1);
        % 
        radial2d   = @(u,sigmax) (X-u(1)).^2/(2*sigmax^2) + (Y-u(2)).^2/(2*sigmax^2);
        gaussian2d = @(u,sigmax)exp( -radial2d(u,sigmax) );
        T = gaussian2d( direction(rho,theta),sigmax) ...
            + gaussian2d(-direction(rho,theta),sigmax);
    case 'gabor-angle'
        rho = getoptions(options, 'rho', 0, 1);
        srho = getoptions(options, 'srho', 0, 1);
        theta = getoptions(options, 'theta', 0, 1);
        stheta = getoptions(options, 'stheta', 0, 1);
        %
        R = sqrt(X.^2+Y.^2);
        Theta = atan2( Y,X );
        % Angular control
        circ = @(u)min(min(abs(u),abs(u+2*pi)),abs(u-2*pi));
        Angular = @(theta,stheta) ...
              gauss(circ(Theta-theta),stheta) ...
            + gauss(circ(Theta-theta+pi),stheta);
        % 
        Radial = @(rho,srho)gauss( R-rho, srho );      
        %
        T = Radial(rho,srho) .* Angular(theta,stheta);
    case 'matern'
        theta = getoptions(options, 'theta', 0, 1); 
        eta = getoptions(options, 'eta', 0, 1);
        alpha = getoptions(options, 'alpha', 0, 1);
        epsilon = getoptions(options, 'epsilon', 1);
        rot  = @(theta)[cos(theta) sin(theta); -sin(theta) cos(theta)];
        tens = @(theta,eta)rot(theta)*diag([1 eta])*rot(-theta);
        TensApply = @(T)T(1,1)*X.^2 + T(2,2)*Y.^2 + 2*T(1,2)*X.*Y;
        %
        T = ( epsilon^2 + TensApply(tens(theta,eta)) ).^(-alpha/2);
    otherwise
        error('Unknown mode');
end
        
T = equalize(T);
