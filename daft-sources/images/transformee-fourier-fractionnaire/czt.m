function g = czt(x, k, w, a)
%CZT	Chirp z-transform.
%	G = CZT(X, M, W, A) returns the M-element z-transform of data X,
%	where M, W and A are scalars which specify the contour in the z-plane
%	on which the z-transform is computed.  M is the length of the transform,
%	W is the complex ratio between points on the contour, and A is the
%	complex starting point.  More explicitly, the contour in the z-plane
%	(a spiral or "chirp" contour) is described by
%	    z = A * W.^(-(0:M-1))
%
%	The parameters M, W, and A are optional; their default values are 
%	M = length(X), W = exp(-j*2*pi/M), and A = 1.  These defaults
%	cause CZT to return the z-transform of X at equally spaced points
%	around the unit circle, equivalent to FFT(X).
%
%	If X is a matrix, the chirp z-transform operation is applied to each
%	column.
%
%	See also FFT, FREQZ.

%	Author(s): C. Denham, 1990.
%		   J. McClellan, 7-25-90, revised
%		   C. Denham, 8-15-90, revised
%		   T. Krauss, 2-16-93, updated help
%	Copyright (c) 1984-94 by The MathWorks, Inc.
%       $Revision: 1.1 $  $Date: 2002/11/05 01:18:42 $

%	References:
%	  [1] Oppenheim, A.V. & R.W. Schafer, Discrete-Time Signal
%	      Processing,  Prentice-Hall, pp. 623-628, 1989.
%	  [2] Rabiner, L.R. and B. Gold, Theory and Application of
%	      Digital Signal Processing, Prentice-Hall, Englewood
%	      Cliffs, New Jersey, pp. 393-399, 1975.

[m, n] = size(x); oldm = m;
if m == 1, x = x(:); [m, n] = size(x); end

if nargin < 2, k = length(x); end
if nargin < 3, w = exp(-sqrt(-1) .* 2 .* pi ./ k); end
if nargin < 4, a = 1; end

if any([size(k) size(w) size(a)]~=1),
	error('Inputs M, W and A must be scalars.')
end

%------- Length for power-of-two fft.

nfft = 1;
while nfft < (m + k - 1), nfft = 2 .* nfft; end

%------- Premultiply data.

kk = ( (-m+1):(k-1) ).';
kk2 = (kk .^ 2) ./ 2;
ww = w .^ (kk2);   % <----- Chirp filter is 1./ww
aa = a .^ (  -kk( m:(2*m-1) )   );
y = x .* [   ( aa.*ww( m:(2*m-1) )  ) * ones(1, n)   ];

%------- Fast convolution via FFT.

fy = fft(  y, nfft );
fv = fft( 1 ./ ww, nfft );   % <----- Chirp filter.
fy = fy .* ( fv * ones(1, n) );
g  = ifft( fy );

%------- Final multiply.

g = g( m:(m+k-1), : ) .* [ ww( m:(m+k-1) )  *  ones(1, n)   ];

if oldm == 1, g = g.'; end
