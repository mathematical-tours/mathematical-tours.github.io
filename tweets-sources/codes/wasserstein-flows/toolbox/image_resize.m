function M1 = image_resize(M,p1,q1,r1)

% image_resize - resize an image using bicubic interpolation
%
%   M1 = image_resize(M,nx,ny,nz);
% or
%   M1 = image_resize(M,newsize);
%
%   Works for 2D, 2D 2 or 3 channels, 3D images.
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin==2
    % size specified as an array
    q1 = p1(2);
    if length(p1)>2
        r1 = p1(3);
    else
        r1 = size(M,3);
    end
    p1 = p1(1);
end

if nargin<4
    r1 = size(M,3);
end


if ndims(M)>=3
    % RVB image
    M1 = zeros(p1,q1, size(M,3), size(M,4));
    for i=1:size(M,3)
        for j=1:size(M,4)
            M1(:,:,i,j) = image_resize(M(:,:,i,j), p1, q1);
        end
    end
    return;
% elseif ndims(M)==3
%     p = size(M,1);
%     q = size(M,2);
%     r = size(M,3);
%     [Y,X,Z] = meshgrid( (0:q-1)/(q-1), (0:p-1)/(p-1), (0:r-1)/(r-1)  );
%     [YI,XI,ZI] = meshgrid( (0:q1-1)/(q1-1), (0:p1-1)/(p1-1), (0:r1-1)/(r1-1) );
%     M1 = interp3( Y,X,Z, M, YI,XI,ZI ,'cubic');
%     return;
end

p = size(M,1);
q = size(M,2);
[Y,X] = meshgrid( (0:q-1)/(q-1), (0:p-1)/(p-1) );
[YI,XI] = meshgrid( (0:q1-1)/(q1-1), (0:p1-1)/(p1-1) );
M1 = interp2( Y,X, M, YI,XI ,'cubic');
