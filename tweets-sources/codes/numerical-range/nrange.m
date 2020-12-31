function [rerange, imrange, nradius] = nrange(A, P, C, F);
%   NRANGE Numerical range of a square matrix
%   [RERANGE, IMRANGE, NRADIUS] = NRANGE(A, P) requires a data square
%   matrix A and the number of points P to be plotted on the boundary
%   of the numerical range of A.
%   It returns the real part and the imaginary part of points of the 
%   boundary of the numerical range of A, together with the numerical
%   radius. It also draws the numerical range.
%
%   [RERANGE, IMRANGE, NRADIUS] = NRANGE(A, P, C) fills the drawing
%   of the numerical range with the specified color C. C can be any of
%   the following: 'r', 'm', 'y', 'g', 'c', 'b', 'k'. The default is 'b'.
%
%   [RERANGE, IMRANGE, NRADIUS] = NRANGE(A, P, C, F) indicates if fill
%   is desired or not. F can be 'y' (yes) or 'n' (no). The default is 'y'.

% Reference: The algorithm used is described in the following book.
% Gustafson K.E. and Rao D.K.M., Numerical Range, Springer-Verlag,
% New York, 1997. See page 137.

% Copyright:
% Dora Matache and Valentin Matache
% University of Nebraska Omaha
% Department of Mathematics
% Year 2001

if nargin == 1
    P = 500;
    C = 'b';
    F = 'y';
end

if nargin == 2
    C = 'b';
    F = 'y';
end

if nargin == 3
    F = 'y';
end

[m, n] = size(A);
if m ~= n
    error('Requires a square matrix.');
end

count = 1;
for theta=0:P
    T=(cos(theta*2*pi/P)+i*sin(theta*2*pi/P))*A;
    ReT=.5*(T+T');
    lambda=max((eig(ReT)));
    [eigvec,eigval]=eig(ReT);
    for j=1:n
        value=eigval(j,j);
        if value>=lambda;
            u=eigvec(:,j)/norm(eigvec(:,j));
            z=u'*A*u;
            x(count)=real(z);
            y(count)=imag(z);
            count = count + 1;
        end
    end
end
rerange = x(1:count-1);
imrange = y(1:count-1);
% disp('The numerical radius is:');
nradius = max(sqrt(x.^2+y.^2));
if F == 'n'
    if A == A(1,1)*eye(n)
        plot(real(A(1,1)),imag(A(1,1)),'color', C);
    else
        plot(rerange, imrange, 'color', C);
        axis equal
    end
else
    if A == A(1,1)*eye(n)
        plot(real(A(1,1)),imag(A(1,1)),'color', C);
    else
        fill(rerange, imrange, .5 + .5*C, 'LineWidth', 2, 'EdgeColor', C);
        axis equal
    end
end