function [II,DA] = inv_2x2(AA)
%INV_2X2 calc. the inverses for a block of 2-by-2 matrices.
%   [IA,DA] = INV_2X2(AA) returns a set of 'inverses' IA and
%   an array of determinants DA for the set of 2-by-2 linear 
%   systems in AA. SIZE(AA), SIZE(IA) = [2,2,N], where N is 
%   the number of linear systems. DA is an N-by-1 array of 
%   determinant values. Note that each IA(:,:,K) is an 'inc-
%   omplete inverse DET(A(:,:,K)) * A(:,:,K)^(-1) to improve
%   numerical robustness. To solve a linear system, A*X = B,
%   compute (I*B)./D, given D is non-zero. 
%
%   See also INV_3X3

%   Darren Engwirda : 2018 --
%   Email           : de2363@columbia.edu
%   Last updated    : 03/05/2018

%---------------------------------------------- basic checks    
    if (  ~isnumeric(AA))
        error('inv_2x2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(AA) ~= +3 )
        error('inv_2x2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(AA,1)~= +2 || ...
        size(AA,2)~= +2 )
        error('inv_2x2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%---------------------------------------------- build inv(A)
    II = zeros(size (AA)) ;

    DA = det_2x2(AA) ;
    
    II(1,1,:) = AA(2,2,:) ;
    II(2,2,:) = AA(1,1,:) ;
    II(1,2,:) =-AA(1,2,:) ;
    II(2,1,:) =-AA(2,1,:) ;
    
end

function [DA]    = det_2x2(AA)

    DA = ...
    AA(1,1,:).* AA(2,2,:) - ...
	AA(1,2,:).* AA(2,1,:) ;

end



