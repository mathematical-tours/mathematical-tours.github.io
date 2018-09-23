function [bb] = pwrbal2(pp,pw,tt)
%PWRBAL2 compute the ortho-balls associated with a 2-simplex
%triangulation embedded in R^2 or R^3.
%   [BB] = PWRBAL2(PP,PW,TT) returns the set of power balls
%   associated with the triangles in [PP,TT], such that BB = 
%   [XC,YC,RC.^2]. PW is a vector of vertex weights.

%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 02/05/2018

%---------------------------------------------- basic checks    
    if ( ~isnumeric(pp) || ...
         ~isnumeric(pw) || ...
         ~isnumeric(tt) )
        error('pwrbal2:incorrectInputClass' , ...
            'Incorrect input class.');
    end
    
%---------------------------------------------- basic checks
    if (ndims(pp) ~= +2 || ...
        ndims(pw) ~= +2 || ...
        ndims(tt) ~= +2 )
        error('pwrbal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
    if (size(pp,2) < +2 || ...
            size(pp,1)~= size(pw,1) || ...
                size(tt,2) < +3 )
        error('pwrbal2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end

    switch (size(pp,2))
        case +2
    %-------------------------------------------- alloc work
        AA = zeros(2,2,size(tt,1)) ;
        Rv = zeros(2,1,size(tt,1)) ;
        bb = zeros(size(tt,1),3,1) ;
        
    %-------------------------------------------- lhs matrix     
        ab = pp(tt(:,2),:)-pp(tt(:,1),:);
        ac = pp(tt(:,3),:)-pp(tt(:,1),:);
        
        AA(1,1,:) = ab(:,1) * +2.0 ;
        AA(1,2,:) = ab(:,2) * +2.0 ;      
        AA(2,1,:) = ac(:,1) * +2.0 ;
        AA(2,2,:) = ac(:,2) * +2.0 ;
       
    %-------------------------------------------- rhs vector    
        Rv(1,1,:) = sum(ab.*ab,2) ...
        - ( pw(tt(:,2)) - pw(tt(:,1)) ) ;
        
        Rv(2,1,:) = sum(ac.*ac,2) ...
        - ( pw(tt(:,3)) - pw(tt(:,1)) ) ;
        
    %-------------------------------------------- solve sys.
       [II,dd] = inv_2x2(AA) ;
        
        bb(:,1) = ( ...
            II(1,1,:).*Rv(1,1,:) + ...
            II(1,2,:).*Rv(2,1,:) ) ...
            ./ dd ;
                
        bb(:,2) = ( ...
            II(2,1,:).*Rv(1,1,:) + ...
            II(2,2,:).*Rv(2,1,:) ) ...
            ./ dd ;
          
        bb(:,1:2) = ...
            pp(tt(:,1),:) + bb(:,1:2) ;
            
    %-------------------------------------------- mean radii
        r1 = sum( ...
        (bb(:,1:2)-pp(tt(:,1),:)).^2,2) ;
        r2 = sum( ...
        (bb(:,1:2)-pp(tt(:,2),:)).^2,2) ;
        r3 = sum( ...
        (bb(:,1:2)-pp(tt(:,3),:)).^2,2) ;

        r1 = r1 - pw(tt(:,1));
        r2 = r2 - pw(tt(:,2));
        r3 = r3 - pw(tt(:,3));

        bb(:,3) = ( r1+r2+r3 ) / +3.0 ;

        case +3
    %-------------------------------------------- alloc work
        AA = zeros(3,3,size(tt,1)) ;
        Rv = zeros(3,1,size(tt,1)) ;
        bb = zeros(size(tt,1),4,1) ;
        
    %-------------------------------------------- lhs matrix     
        ab = pp(tt(:,2),:)-pp(tt(:,1),:);
        ac = pp(tt(:,3),:)-pp(tt(:,1),:);
        
        AA(1,1,:) = ab(:,1) * +2.0 ;
        AA(1,2,:) = ab(:,2) * +2.0 ;
        AA(1,3,:) = ab(:,3) * +2.0 ;      
        AA(2,1,:) = ac(:,1) * +2.0 ;
        AA(2,2,:) = ac(:,2) * +2.0 ;
        AA(2,3,:) = ac(:,3) * +2.0 ;
        
        nv = cross(ab,ac) ;
        
        AA(3,1,:) = nv(:,1) * +1.0 ;
        AA(3,2,:) = nv(:,2) * +1.0 ;
        AA(3,3,:) = nv(:,3) * +1.0 ;
        
    %-------------------------------------------- rhs vector    
        Rv(1,1,:) = sum(ab.*ab,2) ...
        - ( pw(tt(:,2)) - pw(tt(:,1)) ) ;
        
        Rv(2,1,:) = sum(ac.*ac,2) ...
        - ( pw(tt(:,3)) - pw(tt(:,1)) ) ;
        
    %-------------------------------------------- solve sys.
       [II,dd] = inv_3x3(AA) ;
        
        bb(:,1) = ( ...
            II(1,1,:).*Rv(1,1,:) + ...
            II(1,2,:).*Rv(2,1,:) + ...
            II(1,3,:).*Rv(3,1,:) ) ...
            ./ dd ;
                
        bb(:,2) = ( ...
            II(2,1,:).*Rv(1,1,:) + ...
            II(2,2,:).*Rv(2,1,:) + ...
            II(2,3,:).*Rv(3,1,:) ) ...
            ./ dd ;
            
        bb(:,3) = ( ...
            II(3,1,:).*Rv(1,1,:) + ...
            II(3,2,:).*Rv(2,1,:) + ...
            II(3,3,:).*Rv(3,1,:) ) ...
            ./ dd ;
          
        bb(:,1:3) = ...
            pp(tt(:,1),:) + bb(:,1:3) ;
                
    %-------------------------------------------- mean radii
        r1 = sum( ...
        (bb(:,1:3)-pp(tt(:,1),:)).^2,2) ;
        r2 = sum( ...
        (bb(:,1:3)-pp(tt(:,2),:)).^2,2) ;
        r3 = sum( ...
        (bb(:,1:3)-pp(tt(:,3),:)).^2,2) ;

        r1 = r1 - pw(tt(:,1));
        r2 = r2 - pw(tt(:,2));
        r3 = r3 - pw(tt(:,3));

        bb(:,4) = ( r1+r2+r3 ) / +3.0 ;
        
    otherwise

    error('pwrbal2:unsupportedDimension' , ...
            'Dimension not supported.') ;
    
    end

end



