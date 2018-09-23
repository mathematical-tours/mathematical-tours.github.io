function D = pointCloud2TDF(V,X,Y,Z)
% Given a point cloud, compute a voxel grid of TDF values (warning: slow!)
% Works fine for small point clouds - uses in-house Matlab functions
%
% ---------------------------------------------------------
% Copyright (c) 2016, Andy Zeng
% 
% This file is part of the 3DMatch Toolbox and is available 
% under the terms of the Simplified BSD License provided in 
% LICENSE. Please retain this notice and LICENSE if you use 
% this file (or any portion of it) in your project.
% ---------------------------------------------------------

% Compute voxel grid
% [X,Y,Z] = ndgrid((xRange(1)+voxelSize/2):voxelSize:(xRange(2)-voxelSize/2), ...
%                             (yRange(1)+voxelSize/2):voxelSize:(yRange(2)-voxelSize/2), ...
%                             (zRange(1)+voxelSize/2):voxelSize:(zRange(2)-voxelSize/2));

% Build KD-tree and do 1-NN search                         
modelKDT = KDTreeSearcher(V');
[I,D] = knnsearch(modelKDT,[X(:),Y(:),Z(:)]);

% Reshape values into voxel grid
D = reshape(D,size(X));

% Apply truncation
% D = 1.0 - min(D./(voxelSize*voxelMargin),1.0);

end