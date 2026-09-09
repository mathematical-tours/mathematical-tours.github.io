function [vertex,face] = read_off(filename)

% read_off - read data from OFF file.
%
%   Copyright (c) 2003 Gabriel Peyré

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end
buffer = fgetl(fid);
if( buffer~='OFF' )
    error('The file is not a valid OFF one.');
    return;
end
buffer = fgetl(fid);
[nvert,buffer] = strtok(buffer);
nvert = str2num(nvert);
[nface,buffer] = strtok(buffer);
nface = str2num(nface);

vertex = zeros(nvert,3);
face = zeros(nface,3);

str = ['Load mesh ', filename, ' : ', num2str(nvert), ' vertices, ', num2str(nface), ' faces.'];
disp(str);
disp('Loading vertices.');
str = strrep(str,'_','\_'); % problem when displaying '_' (LaTeX special char)
str = strrep(str,'\','/'); % problem when displaying '\' (LaTeX special char)
h = waitbar(0,str);
% read vertex
for i=1:nvert
    waitbar(0.5*i/nvert);
    buffer = fgetl(fid);
    [x,buffer] = strtok(buffer);
    [y,buffer] = strtok(buffer);
    [z,buffer] = strtok(buffer);
    vertex(i,:) = [str2double(x) str2double(y) str2double(z)];
end
% read face
disp('Loading faces.');
for i=1:nface
    waitbar(0.5*i/nface+0.5);    
    buffer = fgetl(fid);
    [n,buffer] = strtok(buffer);
    if( str2num(n)~=3 )
        error('Mesh is not a valid triangulation.');
        return;
    end;
    [a,buffer] = strtok(buffer);
    [b,buffer] = strtok(buffer);
    [c,buffer] = strtok(buffer);
    face(i,:) = [str2num(a) str2num(b) str2num(c)];
end
face = face+1;

close(h);