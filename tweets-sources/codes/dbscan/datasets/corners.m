function data = corners(N, scale, gapwidth, cornerwidth)

    if nargin < 1
        N = 1000;
    end
    if mod(N,8) ~= 0
        N = round(N/8) * 8;
    end

    if nargin < 2
        scale = 10;
    end
    if nargin < 3
        gapwidth = 2;
    end   
    if nargin < 4
        cornerwidth = 2;
    end

    perCorner = N/4;

    xplusmin = [ones(perCorner,1); -1*ones(perCorner,1); ones(perCorner,1); -1*ones(perCorner,1)];
    yplusmin = [ones(perCorner,1); -1*ones(2*perCorner,1); ones(perCorner,1)];
    
    horizontal = [xplusmin(1:2:end) * gapwidth + xplusmin(1:2:end) * scale .* rand(N/2,1), ...
                  yplusmin(1:2:end) * gapwidth + cornerwidth * yplusmin(1:2:end) .* rand(N/2,1), ...
                  floor((0:N/2-1)'/(perCorner*.5))];
       
    vertical = [xplusmin(2:2:end) * gapwidth + cornerwidth * xplusmin(2:2:end) .* rand(N/2,1), ...
                yplusmin(2:2:end) * gapwidth + yplusmin(2:2:end) * scale .* rand(N/2,1), ...
                floor((0:N/2-1)'/(perCorner*.5))];
    
    data=  [horizontal; vertical];

end