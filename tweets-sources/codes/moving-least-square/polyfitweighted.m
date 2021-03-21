function p = polyfitweighted(x,y,n,w)
% polyfitweighted.m 
% -----------------
%
% Find a least-squares fit of 1D data y(x) with an nth order 
% polynomial, weighted by w(x).
%
% By S.S. Rogers (2006), based on polyfit.m by The MathWorks, Inc. - see doc
% polyfit for more details.
%
% Usage
% -----
%
% P = polyfitweighted(X,Y,N,W) finds the coefficients of a polynomial 
% P(X) of degree N that fits the data Y best in a least-squares sense. P 
% is a row vector of length N+1 containing the polynomial coefficients in
% descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1). W is
% a vector of weights. 
%
% Vectors X,Y,W must be the same length.
%
% Class support for inputs X,Y,W:
%    float: double, single
%

% The regression problem is formulated in matrix format as:
%
%    yw = V*p    or
%
%          3    2
%    yw = [x w  x w  xw  w] [p3
%                            p2
%                            p1
%                            p0]
%
% where the vector p contains the coefficients to be found.  For a
% 7th order polynomial, matrix V would be:
%
% V = [w.*x.^7 w.*x.^6 w.*x.^5 w.*x.^4 w.*x.^3 w.*x.^2 w.*x w];

if ~isequal(size(x),size(y),size(w))
    error('X and Y vectors must be the same size.')
end

x = x(:);
y = y(:);
w = w(:);


% Construct weighted Vandermonde matrix.
%V(:,n+1) = ones(length(x),1,class(x));
V(:,n+1) = w;
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
ws = warning('off','all'); 
p = R\(Q'*(w.*y));    % Same as p = V\(w.*y);
warning(ws);
if size(R,2) > size(R,1)
   warning('polyfitweighted:PolyNotUnique', ...
       'Polynomial is not unique; degree >= number of data points.')
elseif condest(R) > 1.0e10
    if nargout > 2
        warning('polyfitweighted:RepeatedPoints', ...
            'Polynomial is badly conditioned. Remove repeated data points.')
    else
        warning('polyfitweighted:RepeatedPointsOrRescale', ...
            ['Polynomial is badly conditioned. Remove repeated data points\n' ...
            '         or try centering and scaling as described in HELP POLYFIT.'])
    end
end
p = p.';          % Polynomial coefficients are row vectors by convention.
