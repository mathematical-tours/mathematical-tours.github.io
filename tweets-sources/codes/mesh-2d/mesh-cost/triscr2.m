function [tscr] = triscr2(pp,tt)
%TRISCR2 calc. area-len. ratios for triangles in a 2-simplex
%triangulation in the two-dimensional plane.
%   [SCR2] = TRISCR2(VERT,TRIA) returns the area-len. ratios
%   where SCR2 is a T-by-1 vector, VERT is a V-by-2 array of 
%   XY coordinates, and TRIA is a T-by-3 array of
%   vertex indexing, where each row defines a triangle, such 
%   that VERT(TRIA(II,1),:), VERT(TRIA(II,2),:) and VERT(
%   TRIA(II,3),:) are the coordinates of the II-TH triangle.
%
%   See also TRIAREA, TRIANG2, TRIBAL2

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 17/01/2017

%--------------------------- compute signed area-len. ratios
    scal = 4.0 * sqrt(3.0) / 3.0;

    area = triarea(pp,tt) ;             % also error checks!

    lrms = sum((pp(tt(:,2),:)        ...
              - pp(tt(:,1),:)).^2,2) ...
         + sum((pp(tt(:,3),:)        ...
              - pp(tt(:,2),:)).^2,2) ...
         + sum((pp(tt(:,3),:)        ...
              - pp(tt(:,1),:)).^2,2) ;

    lrms =(lrms / 3.0) .^ 1.00 ;

    tscr = scal * area ./ lrms ;

end



