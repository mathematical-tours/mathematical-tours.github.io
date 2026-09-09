function xy_spec1 = rectify_embedding(xy,xy_spec)

% rectify_embedding - try to match the embeding
%   with another one.
%   Use the Y coord to select 2 base points.
%
%   xy_spec = rectify_embedding(xy,xy_spec);
%
%   Copyright (c) 2003 Gabriel Peyré

I = find( xy(:,2)==max(xy(:,2)) );
n1 = I(1); 
I = find( xy(:,2)==min(xy(:,2)) );
n2 = I(1);
v1 = xy(n1,:)-xy(n2,:);
v2 = xy_spec(n1,:)-xy_spec(n2,:);
theta = acos( dot(v1,v2)/sqrt(dot(v1,v1)*dot(v2,v2)) );
theta = theta * sign( det([v1;v2]) );
M = [cos(theta) sin(theta); -sin(theta) cos(theta)];
xy_spec1 = (M*xy_spec')';