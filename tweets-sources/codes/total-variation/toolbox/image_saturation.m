function u = image_saturation(f1,rho)

% image_saturation - saturate an image in [-rho,1+rho] and then rescale to [0,1]
%
%   u = image_saturation(f1,rho);
%
% Copyright (c) 2016 Gabriel Peyre

u = clamp(f1,-rho,1+rho); 
%
i = find(u==min(u(:)));
u(i(1))=-rho; 
%
i = find(u==max(u(:)));
u(i(1))=rho; 
u = rescale(u);

end
