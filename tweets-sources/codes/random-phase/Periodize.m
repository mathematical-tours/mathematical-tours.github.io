function p = Periodize(f)

n = size(f,1);

z = zeros(1,1,size(f,3));
fe = [z, f(1,:,:), z; f(:,1,:), f, f(:,end,:); z, f(end,:,:), z];

%%
% Compute the inner-Laplacian \(d = \Delta_i f\) as the usual Laplacian of the
% extended image \(\Delta f_e\).

laplacian = @(x)4*x - ( circshift(x,[0 1]) + circshift(x,[1 0]) + circshift(x,[-1 0]) + circshift(x,[0 -1]) );
d = laplacian(fe);
d = d(2:end-1, 2:end-1, :); 

[X Y] = meshgrid(0:n-1, 0:n-1);
U = 4 - 2*cos(2.*X*pi/n) - 2*cos(2.*Y*pi/n);

%%
% Inverse the Laplacian.

P = fft2(d)./repmat(U, [1 1 size(f,3)]);
P(1,1,:) = sum(sum(f,1),2);
p = real(ifft2(P));

end
