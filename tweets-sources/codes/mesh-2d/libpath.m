function libpath
%LIBPATH a helper function to set-up MATLAB's path statement
%for MESH2D. 
%
%   See also REFINE2, SMOOTH2, TRIDEMO

%-----------------------------------------------------------
%   Darren Engwirda : 2017 --
%   Email           : de2363@columbia.edu
%   Last updated    : 08/07/2018
%-----------------------------------------------------------
%

%------------------------------------ push path to utilities
    
    filename = mfilename('fullpath') ;
    filepath = fileparts( filename ) ;
    
    addpath([filepath,'/aabb-tree']) ;
    addpath([filepath,'/geom-util']) ;
    addpath([filepath,'/hfun-util']) ;
    addpath([filepath,'/hjac-util']) ;
    addpath([filepath,'/mesh-ball']) ;
    addpath([filepath,'/mesh-cost']) ;
    addpath([filepath,'/mesh-file']) ;
    addpath([filepath,'/mesh-util']) ;
    addpath([filepath,'/poly-test']) ;

end



