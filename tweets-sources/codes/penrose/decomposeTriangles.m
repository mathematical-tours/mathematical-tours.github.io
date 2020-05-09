function out = decomposeTriangles(in)
%decomposeTriangles Decompose a table of A, A', B, and B' triangles.
%
%   decomposeTriangles decomposes each triangle in the input table and
%   returns all the decomposed triangles together as a new table.
%
%   EXAMPLE
%   Starting with an A triangle, successively apply triangle
%   decomposition three times and display the result.
%
%       t = aTriangle(1i,0,[]);
%       for k = 1:3
%           t = decomposeTriangles(t);
%       end
%       showLabeledTriangles(t)

%   Copyright 2018 The MathWorks, Inc.

c = cell(height(in),1);

for k = 1:height(in)
    t_k = in(k,:);
    
    switch in.Type(k)
        case "A"
            c{k} = decomposeATriangle(t_k);
            
        case "Ap"
            c{k} = decomposeApTriangle(t_k);
            
        case "B"
            c{k} = decomposeBTriangle(t_k);
            
        case "Bp"
            c{k} = decomposeBpTriangle(t_k);
    end
end

out = cat(1,c{:});
