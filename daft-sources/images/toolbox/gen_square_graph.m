function [A,xy] = gen_square_graph( n, connectivity )

% gen_square_graph - génère un graphe d'une grille carrée, i.e.
%   de (Z/nZ)^2
%   associé à la partie génératrice S = { [-1,0], [1,0], [0,1], [0,-1] } 
%   pour 'connectivity==4' et S = { [-1,-1], [1,1], [-1,1], [1,-1], [-1,0], [1,0], [0,1], [0,-1] } 
%
%   A : matrice d'adjacence de taille n x n
%   xy : matrice de position de taille n x 2.
%
%   Copyright (c) 2003 Gabriel Peyré

if nargin<2
    connectivity = 4;
end

A = zeros(n^2,n^2);
xy = zeros(n^2,2);
h = 1/(n-1);

for i=0:n-1
for j=0:n-1
    k = i+n*j;
    xy(k+1,1) = i*h;
    xy(k+1,2) = j*h;
    if i<n-1
        A( k+1, k+2 ) = 1;
    end
    if i>0
        A( k+1, k ) = 1;
    end
    if j<n-1
        A( k+1, k+1+n ) = 1;
    end
    if j>0
        A( k+1, k+1-n ) = 1;
    end
    if connectivity==8
        if (i<n-1) & (j<n-1)
            A( k+1, k+2+n ) = 1;
        end
        if (i<n-1) & (j>0)
            A( k+1, k+2-n ) = 1;
        end
        if (i>0) & (j>0)
            A( k+1, k-n ) = 1;
        end
        if (i>0) & (j<n-1)
            A( k+1, k+n ) = 1;
        end
    end
end
end