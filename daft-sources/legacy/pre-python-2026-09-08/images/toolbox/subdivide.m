function [vertex1,face1] = subdivide(vertex,face,nsub)

% subdivide - perform a 1:4 subdivision.
%   Subdivide each triangle into 4 smaller triangles.
%   
%   [vertex1,face1] = subdivide(vertex,face,nsub);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    error('Not enough arguments');
end
if nargin==2
    nsub=1;
end

if nsub>1
    % special case for multi-subdivision
    vertex1 = vertex;
    face1 = face;
    for i = 1:nsub
         [vertex1,face1] = subdivide(vertex1,face1,1);
    end
    return;    
end

nface = size(face,1);
nvert = max(max(face));

% this is a hesh-table-like structure to store new vertices
new_verts.test = 0;

m = zeros(3,1); % index of the mid vertices
face1 = [];
for i=1:nface
    f = face(i,:);
    % create new vertices
    for k=1:3
        v1 = f(k);
        v2 = f(mod(k,3)+1);
        str = sprintf('I%dI%dI',min(v1,v2),max(v1,v2)); % a hash string for this edge
        if ~isfield(new_verts,str)
            % we must create a new vertex
            p = ( vertex(v1,:)+vertex(v2,:) )/2;    % position of new vertex
            vertex = [vertex; p];
            m(k) = size(vertex,1);
            new_verts = setfield(new_verts,str,m(k));  % assign to hash table
        else
            m(k) = getfield(new_verts,str);
        end
    end
    % create new face
    face1 = [face1; [m(1),m(2),m(3)] ];     % central face
    face1 = [face1; [f(1),m(1),m(3)] ];
    face1 = [face1; [m(1),f(2),m(2)] ];
    face1 = [face1; [m(2),f(3),m(3)] ];
end

vertex1 = vertex;