function lap = compute_geometric_laplacian(vertex,face,type,ring)

% compute_laplacian - return the combinatorial laplacian
%   of a given triangulation.
%
%   compute_geometric_laplacian(vertex,face,type,ring);
%
%   Type is either : 
%       - 'combinatorial' : combinatorial laplacian, doesn't take into
%       acount geometry.
%       - 'shape_preserving' : Floater "Shape Preserving" weights, 
%       
%
%       Reference: 
%           M.S.Floater and K.Hormann
%           Recent Advances in Surface Parameterization
%           Multiresolution in geometric modelling
%           <http://vcg.isti.cnr.it/~hormann/papers/survey.pdf>
%
%   Copyright (c) 2003 Gabriel Peyré

nface = length(face);
nvert = max(max(face));
lap = zeros(nvert,nvert);



if strcmp(lower(type),'combinatorial')
    lap = compute_laplacian( triangulation2adjacency(face) );
    
elseif strcmp(lower(type),'conformal') || strcmp(lower(type),'authalic')
    if nargin<4
        disp('--> Computing 1-ring.');
        ring = compute_1_ring( face );
    end
    disp('--> Computing laplacian.');
    for i=1:nvert
        vi = vertex(i,:);
        r = ring{i};
        if r(end)==-1
            % no circularity
            s = length(r)-1;
            r = [r(1), r(1:(end-1)), r(end-1)];
        else
            % circularity
            s = length(r);
            r = [r(end), r, r(1)];
        end
        % circulate on the 1-ring        
        for x = 2:(s+1)        
            j = r(x);            
            if lap(i,j)==0            
                gche = r(x-1);                
                drte = r(x+1);                
                vj = vertex(j,:);                
                v1 = vertex(gche,:);                
                v2 = vertex(drte,:);                
                % we use cot(acos(x))=x/sqrt(1-x^2)
                if strcmp(lower(type),'conformal')
                    d1 = sqrt(dot(vi-v2,vi-v2));                
                    d2 = sqrt(dot(vj-v2,vj-v2));                
                    if d1>eps && d2>eps                
                        z = dot(vi-v2,vj-v2)/( d1*d2 );                    
                        lap(i,j) = lap(i,j) + z/sqrt(1-z^2);                    
                    end                
                    d1 = sqrt(dot(vi-v1,vi-v1));                
                    d2 = sqrt(dot(vj-v1,vj-v1));                
                    if d1>eps && d2>eps                
                        z = dot(vi-v1,vj-v1)/( d1*d2 );                    
                        lap(i,j) = lap(i,j) + z/sqrt(1-z^2);                    
                    end         
                else
                    d1 = sqrt(dot(vi-vj,vi-vj));                
                    d2 = sqrt(dot(v2-vj,v2-vj));                
                    if d1>eps && d2>eps                
                        z = dot(vi-vj,v2-vj)/( d1*d2 );                    
                        lap(i,j) = lap(i,j) + z/sqrt(1-z^2);                    
                    end
                    d1 = sqrt(dot(vi-vj,vi-vj));
                    d2 = sqrt(dot(v1-vj,v1-vj));                 
                    if d1>eps && d2>eps                
                        z = dot(vi-vj,v1-vj)/( d1*d2 );                    
                        lap(i,j) = lap(i,j) + z/sqrt(1-z^2);                    
                    end
                    if d1>eps
                        lap(i,j) = lap(i,j) / (d1*d1);
                    end
                end
                if 0    % uncomment for symmeterization
                if lap(j,i)==0
                    lap(j,i) =  lap(i,j);
                else
                    lap(j,i) =  (lap(j,i)+lap(i,j))/2;
                end
                end
            end            
        end        
    end
            
    for i=1:nvert    
        lap(i,i) = -sum( lap(i,:) );        
    end
else
    error('Unknown type.');        
end
 