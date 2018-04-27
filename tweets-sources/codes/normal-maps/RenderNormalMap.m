function [R,G,Ph] = RenderNormalMap(f,l,v, k, col_g, col_ph, alpha)

% R: rendering
% G: Gouraud term
% Ph: Phong term
% alpha: Gouraud exponent
% k=[k_g,k_ph]: Gouraud and Phong weights
%
% f: normal map
% l: light position
% v: view position

n = size(f,1);
normalize = @(f)f ./ repmat( sqrt(sum(f.^2,3)), [1 1 3] );
rm = @(u) repmat( reshape(u,[1 1 3]), [n n 1] );

% pixel position
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t); Z = zeros(n);
P = cat(3, X,Y,Z);
%
% pixel->view vector
V = normalize( rm(v) - P );
% pixel->light vector
L = normalize( rm(l) - P );
% reflected vector, Q=2*<L,N>N-L
Q = 2*repmat(sum(L.*f,3), [1 1 3]).*f - L;
% Gourau term
G = max(0, sum( L.*f, 3 ) );
% Phong term
Ph = max(0, sum( V.*Q, 3 ) ).^alpha;

% k = k/sum(k);
R = zeros(n,n,3);
for j=1:3
    R(:,:,j) = k(1)*col_g(j)*G+k(2)*col_ph(j)*Ph;
end
R = max(R,0);

end