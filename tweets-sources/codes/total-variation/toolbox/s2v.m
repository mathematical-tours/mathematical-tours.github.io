function v = s2v(s,a)

% s2v - structure array to vector
%
%   v = s2v(s,a);
%
%   v(i) = getfield(s(i), a);
%
%   Copyright (c) 2010 Gabriel Peyre

v = [];
for  i = 1:length(s)
    u = getfield(s(i), a);
    if size(u,1)==1 && size(u,2)==1
        v(i) = u;
    elseif (size(u,1)>1 && size(u,2)==1) || (size(u,2)>1 && size(u,1)==1) 
        v(:,i) = u(:);
    else
        v(:,:,i) = u;
    end
end
