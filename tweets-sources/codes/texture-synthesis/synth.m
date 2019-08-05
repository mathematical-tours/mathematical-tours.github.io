function [Image, Mapping] = synth(rawSample, winsize, newRows, newCols, outpath)

%  Non-parametric Texture Synthesis using Efros & Leung's algorithm  
%  Author: Alex Rubinsteyn (alex.rubinsteyn at gmail)


%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301  USA


% Inputs:
% 'filename': the image file containing the sample image (the texture to grow)
% 'winsize': the edge length of the window to match at each iteration (the window is (winsize x winsize) )
% (newRows, newCols): the size of the output image

% Outputs:
% 'Image': the output image (the synthesized texture)
% 'time': the amount of time it took to perform the synthesis


MaxErrThreshold = 0.1;

% rawSample = im2double(imread(filename)); 
rawSample = im2double(rawSample);

sample =  rawSample;  

[rows, cols, channels] = size(sample); 
windowlessSize = [(rows - winsize + 1) (cols - winsize + 1)];

halfWindow = (winsize - 1) / 2;

npixels = newRows * newCols; 

Image = zeros(newRows, newCols, 3); 
Mapping = zeros(newRows, newCols,3);

red_patches = im2col(sample(:, :, 1), [winsize winsize], 'sliding'); 
green_patches = im2col(sample(:, :, 2), [winsize winsize], 'sliding'); 
blue_patches = im2col(sample(:, :, 3), [winsize winsize], 'sliding'); 


%initialize new texture with a random 3x3 patch from the sample
randRow = ceil(rand() * (rows - 2)); 
randCol = ceil(rand() * (cols - 2));

seedSize = 3; 
seedRows = ceil(newRows/2):ceil(newRows/2)+seedSize-1;
seedCols = ceil(newCols/2):ceil(newCols/2)+seedSize-1;


M = sample(randRow:randRow+seedSize-1, randCol:randCol+seedSize-1, :);
if 0
k = length(seedRows); 
    I = randperm(k*k);
for i=1:3
    m = M(:,:,i); m = reshape(m(I), size(m)); 
    M(:,:,i) = m;
end
end
Image(seedRows, seedCols, :) = M;

Mapping(seedRows, seedCols,1) = (randRow:randRow+seedSize-1)' * ones(1,length(seedCols));
Mapping(seedRows, seedCols,2) = ones(length(seedCols),1) * (randCol:randCol+seedSize-1);

nfilled = seedSize * seedSize; 
filled = repmat(false, [newRows newCols]); 
filled(seedRows, seedCols) = repmat(true, [3 3]); 

gaussMask = fspecial('gaussian',winsize, winsize/6.4);

nskipped = 0; 

while nfilled < npixels    
    progress = false;
    
    [pixelRows, pixelCols] = GetUnfilledNeighbors(filled, winsize);
     
     for i = 1:length(pixelRows)
        pixelRow = pixelRows(i);
        pixelCol = pixelCols(i);
        
        rowRange = pixelRow-halfWindow:pixelRow+halfWindow;
        colRange =  pixelCol - halfWindow:pixelCol + halfWindow;

        deadRows = rowRange < 1 | rowRange > newRows;
        deadCols = colRange < 1 | colRange > newCols; 


        if sum(deadRows) + sum(deadCols) > 0 
            safeRows = rowRange(~deadRows); 
            safeCols = colRange(~deadCols); 

            template = zeros(winsize, winsize, 3); 
            template(~deadRows, ~deadCols, :) = Image(safeRows, safeCols, :); 

            validMask = repmat(false, [winsize winsize]); 
            validMask(~deadRows, ~deadCols) = filled(safeRows, safeCols); 
        else
            template = Image(rowRange, colRange, :);
            validMask = filled(rowRange, colRange); 

        end

       [bestMatches, SSD] = FindMatches(template, validMask, gaussMask, red_patches, green_patches, blue_patches);


        matchIdx = RandomPick(bestMatches);
        matchError = SSD(matchIdx); 

         if matchError < MaxErrThreshold 
             [matchRow, matchCol] = ind2sub(windowlessSize, matchIdx); 
             
             %match coords are at corner of window and need to be offset
             matchRow = matchRow + halfWindow;
             matchCol = matchCol + halfWindow;  

             Mapping(pixelRow, pixelCol,1) = matchRow;
             Mapping(pixelRow, pixelCol,2) = matchCol;
             Image(pixelRow, pixelCol, :) = sample(matchRow, matchCol, :);

             filled(pixelRow, pixelCol) = true;   
             nfilled = nfilled + 1; 
             progress = true;
         else
             nskipped = nskipped + 1; 
         end
    end
    
    progressbar(nfilled,npixels);
    if not(exist('k'))
        k=1;
    end
    I = Image + (1-filled);
    J = Mapping/max(rows, cols) + (1-filled);
    clf; 
    subplot(1,2,1); image(I); axis image; axis off; 
    subplot(1,2,2); image(J); axis image; axis off; 
    drawnow;
    imwrite(rescale(I), [outpath 'anim-' znum2str(k,3) '.png' ]);
    imwrite(rescale(J), [outpath 'map-' znum2str(k,3) '.png' ]);
    k = k+1;
%     disp(sprintf('Pixels filled: %d / %d', nfilled, npixels)); 
%     
%     figure; 
%     subplot(2,1,1); 
%     imshow(filled); 
%     subplot(2,1,2); 
%     imshow(Image); 
    
    if ~progress 
        MaxErrThreshold = MaxErrThreshold * 1.1;
        disp(sprintf('Incrementing error tolerance to %d', MaxErrThreshold)); 
    end
end


%% Get pixels at edge of synthesized texture
 function [pixelRows, pixelCols] = GetUnfilledNeighbors(filled, winsize) 
    border = bwmorph(filled,'dilate')-filled;
    
    [pixelRows, pixelCols] = find(border);
    len = length(pixelRows); 
     
     %randomly permute candidate pixels
     randIdx = randperm(len); 
     pixelRows = pixelRows(randIdx); 
     pixelCols = pixelCols(randIdx); 

     %sort by number of neighbors     
     filledSums = colfilt(filled, [winsize winsize], 'sliding', @sum); 
     numFilledNeighbors = filledSums( sub2ind(size(filled), pixelRows, pixelCols) ); 
     [sorted, sortIndex] = sort(numFilledNeighbors, 1, 'descend');
     
     pixelRows = pixelRows(sortIndex); 
     pixelCols = pixelCols(sortIndex); 

     
%% Pick a random pixel from valid patches
function idx = RandomPick(matches)
    indices = find(matches);
    idx = indices(ceil(rand() * length(indices))); 
 
     
%% Find candidate patches that match template
function [pixelList, SSD] = FindMatches (template, validMask, gaussMask, red_patches, green_patches, blue_patches)
ErrThreshold = 0.3; 

[pixels_per_patch, npatches] = size(red_patches); 

totalWeight = sum(sum(gaussMask(validMask)));

mask = (gaussMask .* validMask) / totalWeight;
mask_vec = mask(:)'; 
 
red = reshape(template(:, :, 1), [pixels_per_patch 1]); 
green = reshape(template(:, :, 2), [pixels_per_patch 1]); 
blue = reshape(template(:, :, 3), [pixels_per_patch 1]);

red_templates = repmat(red, [1 npatches]); 
green_templates = repmat(green, [1 npatches]); 
blue_templates = repmat(blue, [1 npatches]); 

red_dist =  mask_vec * (red_templates - red_patches).^2; 
green_dist = mask_vec * (green_templates - green_patches).^2 ; 
blue_dist = mask_vec * (blue_templates - blue_patches).^2; 

SSD = (red_dist + green_dist + blue_dist); 

pixelList = SSD <= min(SSD) * (1+ErrThreshold);  
