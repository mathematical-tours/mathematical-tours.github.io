function plot_gear_animation(x,y,phi,L, q, gifname)

% plot_gear_animation - display and animate a primal/dual gear pair.
%
%   plot_gear_animation(x,y,phi,L, q, gifname);
%
%   x,y are the two gears
%   phi is the mapping between the angles
%   L is the spacing between the gears
%   q is the number of frame for a complete rotation.
%   gifname is the name of the file to output an animated gif.
%
%   Copyright (c) 2010 Gabriel Peyre


if nargin<5
    q = 100; % number of time step
end
if nargin<6
    gifname = '';
end

n = length(x);
b = max(max(x),max(y));
a1 = max(x); a2 = L+max(y);

t = 2*pi* (0:n-1)'/n;
tx = linspace(0,2*pi,q+1)';
ty = linspace(0,2*pi,n+1)';
ty = interp1( ty, phi, tx);

% initialize for gif display
if not(isempty(gifname))
    fprintf('Beginning gif file capture ... ');
end
im = [];
alreadysaved = 0;

while true
    
    for i=1:q
        clf;
        hold on;
        plot_gear(x, [0 0], 'b',0, t + tx(i));
        plot_gear(y, [L 0], 'r',1, t + ty(i));
        axis([-a1 a2 -b b]);
        set(gcf, 'Color', [1 1 1]);
        drawnow;
        if not(alreadysaved) && not(isempty(gifname))
            f = getframe;
            if i==1
                [im,map] = rgb2ind(f.cdata,256,'nodither');
            else
                im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
            end
        end
    end
    
    if not(alreadysaved) && not(isempty(gifname))
        % save to fle
        imwrite(im,map,gifname,'DelayTime',0,'LoopCount',inf) %g443800
        alreadysaved = 1;
        fprintf('done.\n');
    end

end

