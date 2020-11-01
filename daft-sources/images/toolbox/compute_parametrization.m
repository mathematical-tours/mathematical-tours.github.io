function xy = compute_parametrization(vertex,face,lap_type,boundary_type,ring)

% compute_parametrization - compute a planar parameterization
%   of a given disk-like triangulated manifold.
%
%   compute_parametrization(vertex,face,lap_type,boundary_type,ring);
%
%   'lap_type' is either 'combinatorial', 'conformal' or 'authalic'
%   'boundary_type' is either 'circle', 'square' or 'triangle'

nface = length(face);
nvert = max(max(face));

if nargin<2
    error('Not enough arguments.');
end
if nargin<3
    lap_type = 'conformal';
end
if nargin<4
    boundary_type = 'circle';
end
if nargin<5
    disp('--> Computing 1-ring.');
    ring = compute_1_ring( face );
end

% compute Laplacian
lap = compute_geometric_laplacian(vertex,face,lap_type,ring);

% compute boundary
v = -1;
for i=1:nvert   % first find a starting vertex
    f = ring{i};
    if f(end)<0
        v = i;
        break;
    end
end
if v<0
    error('No boundary found.');
end
boundary = [v];
prev = -1;
while true
    f = ring{v};
    if f(end)>=0
        error('Problem in boundary');
    end
    if f(1)~=prev
        prev = v;
        v = f(1);
    else
        prev = v;
        v = f(end-1);
    end
    if ~isempty( find(boundary==v) )
        % we have reach the begining of the boundary
        if v~=boundary(1)
            warning('Begining and end of boundary doesn''t match.');
        else
            break;
        end
    end
    boundary = [boundary,v];
end

% compute the position of the boundary vertex
nbound = length(boundary);
xy_boundary = zeros(nbound,2);
% compute total length
L = 0;
ii = boundary(end); % last index
for i=boundary
    L = L + norme( vertex(i,:)-vertex(ii,:) );
    ii = i;
end
% compute boundary vertex position
l = 0;
ii = boundary(end); % last index
for k=1:nbound
    i = boundary(k);
    ii = boundary( mod(k,nbound)+1 );   % next point
    t = l/L;
    
    switch lower(boundary_type)
        case 'circle'
            xy_boundary(k,:) = [cos(2*pi*t),sin(2*pi*t)];
        case 'square'
            if t<0.25
                xy_boundary(k,:) = [4*t,0];
            elseif t<0.5
                xy_boundary(k,:) = [1,4*t-1];
            elseif t<0.75
                xy_boundary(k,:) = [3-4*t,1];
            else
                xy_boundary(k,:) = [0,4-4*t];
            end
        case 'triangle'
            if t<1/3
                xy_boundary(k,:) = [3*t,0];
            elseif t<2/3
                xy_boundary(k,:) = (3*t-1)*[0.5,sqrt(2)/2] + (2-3*t)*[1,0];
            else
                xy_boundary(k,:) = (3-3*t)*[0.5,sqrt(2)/2];
            end
        otherwise
            error('Unknown boundary type.');
    end
    
    l = l + norme( vertex(i,:)-vertex(ii,:) );
end

% set up the matrix
M = lap;
for i=boundary
    M(i,:) = zeros(1,nvert);
    M(i,i) = 1;
end
% solve the system
xy = zeros(nvert,2);
for coord = 1:2
    % compute right hand side
    x = zeros(nvert,1);
    x(boundary) = xy_boundary(:,coord);
    xy(:,coord) = M\x;
end

