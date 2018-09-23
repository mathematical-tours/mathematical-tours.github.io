function [same,sloc] = setset2(iset,jset)
%SETSET2 a (much) faster variant of ISMEMBER for edge lists.
%   [IN] = SETSET2(ISET,JSET) returns an I-by-1 array IN,
%   with IN(K) = TRUE if ISET(K,:) is present in JSET. This 
%   routine is essentially an optimised ISMEMBER variant de-
%   signed for processing lists of edge indexing. ISET is an
%   I-by-2 array of "query" edges, JSET is a J-by-2 array of
%   edges to test against.
%
%   See also ISMEMBER 

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 29/01/2017

%---------------------------------------------- basic checks
    if ( ~isnumeric(iset) || ~isnumeric(jset) )
        error('setset2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
  
%---------------------------------------------- basic checks
    if (ndims(iset) ~= +2 || ndims(jset) ~= +2)
        error('setset2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(iset,2)~= +2 || size(jset,2)~= +2)
        error('setset2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

%---------------------------------------------- set v1 <= v2
    iset = sort(iset,2) ;
    jset = sort(jset,2) ;    

%-- this is the slow, but easy-to-undertsand version of what
%-- is happening here...  
  
  % if (nargout == +1)
  % same = ismember(iset,jset,'rows') ;
  % else
  % [same,sloc] = ...
  %        ismember(iset,jset,'rows') ;
    
%-- as above, the 'ROWS' based call to ISMEMBER can be sped
%-- up by casting the edge lists (i.e. pairs of UINT32 valu-
%-- es) to DOUBLE, and performing the sorted queries on vec-
%-- tor inputs!
    if (nargout == +1)
    same       = ismember( ...
        iset*[2^31;1], jset*[2^31;1]) ;
    else
   [same,sloc] = ismember( ...
        iset*[2^31;1], jset*[2^31;1]) ;
    end

end



