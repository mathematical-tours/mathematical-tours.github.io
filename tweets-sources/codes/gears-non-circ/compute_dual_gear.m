function [y,L,phi] = compute_dual_gear(x,K)

% compute_dual_gear - compute the conjugate gear
%
%   [y,L,phi] = compute_dual_gear(x);
%
%   y is the dual to x (radius of the gear as a function of
%       uniformly spaced angles).
%   L is the spacing between the centers of the gears
%   phi is the mapping between the angle, theta_y=phi(theta_x) that should
%       be used to rotate the dual gear.
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<2
    % number of turns
    K = 1;
end

n = length(x);

%%
% Dichotomic search for the separation

evalR= @(L)sum( x ./ (L-x) ) *2*pi/n - 2*pi/K;
niter = 50; % iteration for the search
L1 = max(x)*1.001; L2 = max(x)*5;
if evalR(L1)*evalR(L2)>0
    error('Problem');
end
for i=1:niter
    L = (L1 + L2)/2;
    if evalR(L)<0
        L2 = L;
    else
        L1 = L;
    end
end
L = (L1+L2)/2;

%%
% Change of variables between the anglular cordinates

% forward change of variable
phi = [0; cumsum( x ./ (L-x) ) * 2*pi /n];
phi = phi/phi(end)*2*pi;
% backward change of variable
t = linspace(0,2*pi,n+1)';
phiInv = interp1(phi, t, t);
% new curve
y = L-interp1(t, x([1:end 1]), phiInv);
y = y(1:end-1);

% put phi in correct direction
phi = 2*pi-phi(end:-1:1);

end
