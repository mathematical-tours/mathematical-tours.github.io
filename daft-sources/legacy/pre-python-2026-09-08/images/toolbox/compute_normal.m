function normal = compute_normal(vertex,face)

% compute_normal - compute the normal at each vertex
%   of a triangulation, using mean of the faces normals.
%
%   normal = compute_normal(vertex,face);
%
%   Copyright (c) 2004 Gabriel Peyré

nface = length(face);
nvert = max(max(face));
normal = zeros(nvert,3);

for i=1:nface
    f = face(i,:);
    % compute the normal to the face
    n = cross( vertex(f(2),:)-vertex(f(1),:),vertex(f(3),:)-vertex(f(1),:) );
    n = n/norm(n,'fro');
    for j=1:3
        normal( f(j),: ) = normal( f(j),: ) + n;
    end
end

% normalize
for i=1:nvert
    n = normal(i,:);
    normal(i,:) = n/norm(n,'fro');
end