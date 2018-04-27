function movie_display(f)

% movie_display - display a 3-D array as a movie.
%
%   movie_display(f);
%
%   Copyright (c) 2012 Gabriel Peyre

s = 2.5;
normalize = @(x)rescale( clamp( (x-mean(x(:)))/std(x(:)), -s,s) );
A = normalize(f)*256;

clf;
% stpo button
uicontrol(...
    'Style','pushbutton', 'String', 'Stop',...
    'Units','Normalized', 'Position', [0.4 0.1 0.2 0.1],...
    'Callback', @MyCallback);
% Axes
ax = axes(...
    'Units','Normalized',...
    'OuterPosition', [0 0.2 1 0.8]);


k = 0;
global run;
run = 1;
while run
    k = mod(k,size(A,3))+1;
    image(A(:,:,k)); axis image; axis off;
    colormap gray(256);
    drawnow;
    set(ax, 'ButtonDownFcn', 'get(ax, ''CurrentPoint'')');
end

end

function MyCallback(a,b,c)
global run;
run = 0;
end
