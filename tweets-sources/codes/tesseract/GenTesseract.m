function [C,I] = GenTesseract(d)

% Generate a cube in dimension d (a Tesseract)
% I is the edge connexion matrix.

% 1D cube
C = [0, 1];
I = [1; 2]; 
% 
for k=1:d-1
    n = size(C,2); % should 2^d
    % duplicate the cube in the new dimension
    C = [ ...
        [C;zeros(1,n)], ...
        [C;ones(1,n)], ...
    ];
    % add connection
    J = [1:n; n+1:2*n];
    I = [I, I+n, J];
end
d = size(C,1);
C = 2*C-1;

end