t = bTriangle([],-1,1);
for k = 1:7
    t = decomposeTriangles(t);
end
showTiles(t)
axis([-.25 0.19 0.14 0.48])
axis off

print -dpng -r150 penrose-screen-shot.png

%   Copyright 2018 The MathWorks, Inc.
