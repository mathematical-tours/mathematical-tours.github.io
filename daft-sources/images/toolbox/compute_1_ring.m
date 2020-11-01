function ring = compute_1_ring(face)

% compute_1_ring - compute the 1 ring of each vertex
%   in a triangulation.
%
%   ring = compute_1_ring(face);
%
%   IMPORTANT: each triangle of 'face'
%   must be sorted in increasing order of vertex
%   number.
%
%   Copyright (c) 2004 Gabriel Peyré

nface = size(face,1);
nvert = max(max(face));

% create empty list
for i=1:nvert
    ring{i} = []; 
end

% assert each face is in increasing order
for i=1:nface
    face(i,:) = sort(face(i,:));
end

% first compute the 1 ring not in 
% corect order
for i=1:nface
    f = face(i,:);
    for j=1:3
        k = f(j);
        r = ring{k};
        k1 = f( mod(j,3)+1 );
        k2 = f( mod(j-2,3)+1 );
        if isempty(find(r==k1))
            r = [r k1];
        end
        if isempty(find(r==k2))
            r = [r k2];
        end
        ring{k} = r;
    end
end

h = waitbar(0,'Computing 1-ring');
% sort 1 ring in circular order
for i=1:nvert
    waitbar(i/nvert)
    r = ring{i};
    
    % pop 1st vertex
    gche = r(1);
    drte = r(1);
    verts = r(1);
    r = r(2:length(r));
    
    while gche>=0 || drte>=0
        % find in ring a vertex to complete the triangle <i,gche,?>
        found = 0;
        if gche>0
            for x=r
                f = sort( [i,gche,x] );
                % search in face list
                for j=1:nface
                    if face(j,:)==f
                        gche = x;
                        found = 1;
                        verts = [x,verts];
                        xf = find(r~=x);    % pop vertex
                        r = r(xf);
                        break;
                    end
                end 
                if(found) break; end;
            end
            if(~found) gche=-1; end;
        end
        found = 0;
        if drte>0
            for x=r
                f = sort( [i,drte,x] );
                % search in face list
                for j=1:nface
                    if face(j,:)==f
                        drte = x;
                        found = 1;
                        verts = [verts,x];
                        xf = find(r~=x);    % pop vertex
                        r = r(xf);
                        break;
                    end
                end 
                if(found) break; end;
            end
            if(~found) 
                drte=-1; 
            end;
        end
    end

    if length(verts)~=length(ring{i})
        warning('Problem in 1-ring.');
    end
    
    % test for circularity or not
    f = sort( [i,verts(1),verts(end)] );
    found = 0;
    for j=1:nface
        if face(j,:)==f
            found = 1;
            break;
        end
    end 
    if(~found) 
        verts = [verts,-1]; 
    end;
    
    ring{i} = verts;
end
close(h);