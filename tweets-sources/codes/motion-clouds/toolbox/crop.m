function M = crop(M,n,c)

% crop - crop an image to reduce its size
%
%   M = crop(M,n,c);
%
%   n is the new size of the image
%   c is the center of the grop
%
%   Copyright (c) 2007 Gabriel Peyr?

n0 = size(M);

if nargin<2
    n = round( n0(1:2)/2 );
end
if nargin<3 || isempty(c)
    c = round( n0(1:2)/2 );
end

if isempty(n)
    return;
end

if length(n)==1
    n = [n n];
end
if length(c)==1
    c = [c c];
end

c = round(c);


selx = c(1)-ceil(n(1)/2)+1:c(1)+floor(n(1)/2);
sely = c(2)-ceil(n(2)/2)+1:c(2)+floor(n(2)/2);

selx(selx<1) = []; 
selx(selx>n0(1)) = [];
sely(sely<1) = []; 
sely(sely>n0(2)) = [];

M = M(selx,sely,:);