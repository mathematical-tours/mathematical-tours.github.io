function result = unweightedBilateral(blurImage,weightImage,nBins,imageBlur,intensitySigma)

% Rescale weights to [0,1]
minVal = min(weightImage(:));
maxVal = max(weightImage(:));

% For convenience, expand the range a bit -- this way no pixels equal 0 or 1
range = maxVal-minVal+1e-5;
minVal = minVal - range*.001;
maxVal = maxVal + range*.001;
range = maxVal - minVal;

% Rescale image -- don't worry, we'll undo this!
weightImage = (weightImage - minVal) / range;

% We rescaled the pixels so now we have to rescale the standard deviation
intensitySigma = intensitySigma / range;

% We wish to evaluate f(x) = \int I(y) K1(x,y) K2(I2(x),I2(y)) dy.
% We'll approximate this using Fredo's "signal processing approach."
% Define the following function:
%       g(x,c) = \int I(y) K1(x,y) K2(c,I2(y)) dy
% Notice that f(x) = g(x,I2(x)).  But, for fixed c this is the same as
% applying the blur kernel to the function I(y) K2(c,I2(y))!

% So, we'll sample c values c_j and compute the image g(x,c_j).  Then,
% we'll slice through this array to compute f.

c = linspace(0,1,nBins); % thanks to our rescaling!
slices = cell(nBins,1);

% Construct the function to blur I(y) K2(c_i,I2(y)) and then apply kernel K1
slices = cell(nBins,1);
for i=1:nBins % can be parfor
    unblurred = blurImage .* normpdf(c(i) - weightImage, 0, intensitySigma);
    slices{i} = imageBlur(unblurred);
end

% Use linear interpolation to carry out the slice
result = zeros(size(blurImage));
for i=1:(nBins-1)
    curPixels = (weightImage >= c(i) & weightImage < c(i+1));
    rWeight = (weightImage - c(i)) / (c(i+1)-c(i));
    lWeight = 1-rWeight;
    result(curPixels) = lWeight(curPixels).*slices{i}(curPixels) + ...
                        rWeight(curPixels).*slices{i+1}(curPixels);
end