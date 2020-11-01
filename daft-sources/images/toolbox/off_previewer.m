function off_previewer(file)

% off_previewer - load a file in OFF
%   format and display the model.
%
%   off_previewer(file);
%
%   Use it with no argument if you want
%   to open a file browser.
%
%   Copyright (c) 2004 Gabriel Peyré

h = figure;
if nargin==0
    [f, pathname] = uigetfile({'*.off','OFF Model Files'},'Pick a file');
    file = [pathname,f];
end
[vertex,face] = read_off(file);
cf = 0.8;
ce = 0;
patch('vertices',vertex,'faces',face,'facecolor',[cf cf cf],'edgecolor',[ce ce ce]);
lighting phong;
camlight infinite; 
camproj('perspective');
axis square; axis off; 
title( sprintf('File %s.', f) );
disp('Press any key to exit.');
pause;
close(h);