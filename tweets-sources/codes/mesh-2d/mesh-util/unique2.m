function [set2,imap,jmap] = unique2(set2)
%UNIQUE2 a (much) faster variant of UNIQUE for edge lists.
%   [SET2] = UNIQUE2(SET2) returns the unique elements of
%   the N-by-2 array SET2, equivalent to a call to UNIQUE,
%   with thw 'ROWS' option active.
%   [SET2,IMAP,JMAP] = UNIQUE2(SET2) returns the additional
%   array arguments as per UNIQUE.
%
%   See also UNIQUE 

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 10/07/2018

%---------------------------------------------- basic checks
    if ( ~isnumeric(set2) )
        error('unique2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
  
%---------------------------------------------- basic checks
    if (ndims(set2) ~= +2 || size(set2,2) ~= +2)
        error('unique2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

%------------------------------ unique edges and re-indexing
  %[set2, imap, jmap] = ...
  %     unique(sort(set2, 2), 'rows');
   
%-- as a (much) faster alternative to the 'ROWS' based call
%-- to UNIQUE above, the edge list (i.e. pairs of UINT32 va-
%-- lues) can be cast to DOUBLE, and the sorted comparisons 
%-- performed on vector inputs! 
    if (nargout <=  +2)
    
    set2 = sort(set2,2);
   [stmp,imap     ] = unique(set2*[2^31;1]);  
    set2 = set2(imap,:);
    
    else
    
    set2 = sort(set2,2);
   [stmp,imap,jmap] = unique(set2*[2^31;1]);  
    set2 = set2(imap,:);
    
    end

end



