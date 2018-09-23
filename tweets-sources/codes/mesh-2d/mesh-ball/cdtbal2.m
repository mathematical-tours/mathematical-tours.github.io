function [cc] = cdtbal2(pp,ee,tt)
%CDTBAL2 compute the modified circumballs associated with a 
%constrained 2-simplex Delaunay triangulation in R^2.
%   [CC] = CDTBAL2(PP,EE,TT) returns the smallest enclosing 
%   balls associated with the triangles in [PP,TT], such th-
%   at CC = [XC,YC,RC.^2]. Such balls never lie outside the 
%   boundaries of the associated CDT. See TRICON2 for info-
%   mation regarding the edge array EE.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 01/10/2017

%---------------------------------------------- basic checks    
    if (~isnumeric(pp) || ~isnumeric(ee) || ...
        ~isnumeric(tt) )
        error('cdtbal2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ndims(ee) ~= +2 || ...
        ndims(tt) ~= +2 )
        error('cdtbal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(pp,2)~= +2 || size(ee,2) < +5 || ...
        size(tt,2) < +6 )
        error('cdtbal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

%----------------------------------------- calc. circumballs
    cc = tribal2(pp,tt);
    
%------------------------ replace with face-balls if smaller
    cc = minfac2(cc,pp,ee,tt,1,2,3) ;
    cc = minfac2(cc,pp,ee,tt,2,3,1) ;
    cc = minfac2(cc,pp,ee,tt,3,1,2) ; 
        
end

function [cc] = minfac2(cc,pp,ee,tt,ni,nj,nk)
%MINFAC2 modify the set of circumballs to constrain centres
%to the boundaries of the CDT.
%   [CM] = MINFAC2(CC,PP,EE,TT,NI,NJ,NK) returns the set of
%   modified circmballs CM, where any ball CC lying outside
%   the boundaries of the CDT [PP,EE,TT] is replaced by the
%   edge-centred diametric ball. [NI,NJ] are the local inde-
%   xes associated with an edge to test. NK is the local in-
%   dex of the opposite vertex.

%------------------------------------------------ outer edge          
    EF = ee(tt(:,ni+3),5) > +0 ;

%------------------------------------------------ edge balls
    bc = (pp(tt(EF,ni),:)+pp(tt(EF,nj),:))*.50;
    
%------------------------------------------------ edge radii
    br = sum((bc(:,1:2)-pp(tt(EF,ni),:)).^2,2)...
       + sum((bc(:,1:2)-pp(tt(EF,nj),:)).^2,2);
    br = br * +0.5 ;
              
%------------------------------------------- enclosing radii
    ll = sum((bc(:,1:2)-pp(tt(EF,nk),:)).^2,2);
    
%------------------------------------------- replace if min.
    bi = br >= ll ...
       & br <= cc(EF,3) ;
    ei = find(EF) ; 
    ti = ei  (bi) ;
    
%------------------------------------------- replace is min.
    cc(ti,1:2) = bc(bi,:) ;
    cc(ti,  3) = br(bi,:) ;

end


