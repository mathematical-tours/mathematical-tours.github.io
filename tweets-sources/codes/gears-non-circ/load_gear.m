function x = load_gear(name, n, center, tooth, smoothing)

% load_gear - create default gears
%
%   x = load_gear(name, n, center, tooth, smoothing);
%
%   n is the number of points used for the discretization.
%   center is the coordinate of the center of rotation (detaul is [0 0])
%   tooth gives the parameters for the tooth extrusion.
%
%   Copyright (c) 2010 Gabriel Peyre

theta = (0:n-1)'/n*2*pi;

if nargin<3
    center = [0 0];
end
if nargin<4
    tooth.transition =.2;
    tooth.nbr = 30;
    tooth.height = .05;
end
if nargin<5
    smoothing=0;
end

switch name
    case 'ellipse-focal'
        
        % ellipse, focal
        bmin = 1;
        epsilon = .5;
        x = bmin*(1-epsilon^2)./( 1-epsilon*cos(theta) );
        
    case 'ellipse-centered'
        % ellipse, centered
        bmin = 1;
        bmax = 3;
        x = bmin*bmax./sqrt( (bmin*cos(theta)).^2 + (bmax*sin(theta)).^2 );
        
    case {'random' 'random-strong'}        
        if strcmp(name, 'random')
            bmin = 1; bmax = 1.5;
            p = 6;
        else
            bmin = 1; bmax = 4;
            p = 10;
        end        
        x = zeros(n,1);
        x(1:p) = exp(2i*rand(p,1));
        x = real(ifft(x));
        x = (x-min(x))/(max(x)-min(x));
        x = x*(bmax-bmin)+bmin;
        
    case 'petal'
        x = sin(4*theta)+1.5;        
        
    case 'petal-8'
        x = sin(8*theta).^2+2;
        
    case 'circle'        
        x = theta*0+1;        
        
    case 'cardioid'
        x = 1.2+cos(theta);
        
    case 'engrenage-16'        
        eta = 16;
        x = 10+ abs(cos(theta*eta)).^.2 .* sign(cos(theta*eta));
        
    case 'discont'        
        x = .5*theta/2*pi + 1;
        
        
    case 'discont-2'        
        x = mod(theta/(2*pi),.5)*3 + 1;
        
    case 'square'        
        a = genpolygon_regular(4,n);
        
    case 'triangle'
        a = genpolygon_regular(3,n);        
        
        
    case 'hexagon'        
        a = genpolygon_regular(6,n);
        
end

%% add tooths

if not(exist('a'))
    %% convert to rectangular
    a = gear2cart(x);
end

a(:,1) = a(:,1) + center(1);
a(:,2) = a(:,2) + center(2);

if smoothing>0    
    niter = round(n * smoothing);
    a = smooth(a,niter);
end

if tooth.nbr>0
    q = n*10; % number of point for the profile curve    
    t = mod((0:q-1)'/q, 1/tooth.nbr)*tooth.nbr;
    e = (tooth_profile(t, tooth.transition)) * 2 * tooth.height;
    % normal   
    asmooth = smooth(a,.02);
    u = compute_normal(asmooth);
    s = compute_curvabs(a);
    e1 = interp1( linspace(0,1,q+1), [e;e(1)], s );
    % offset
    a = a + u .* repmat(e1, [1 2]);
end

x = cart2gear(a);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = smooth(a,s)

n = size(a,1);
t = (-n/2:n/2-1)';
h = exp( -(t.^2)/(2*s^2) );
h = h/sum(h);
% recenter the filter for fft use
h1 = fftshift(h);
filter = @(u)real(ifft(fft(h1).*fft(u)));

a(:,1) = filter(a(:,1));
a(:,2) = filter(a(:,2));

return;

if nargin<2
    niter = 1;
end
for i=1:niter
    a = (a + a([2:end 1],:) + a([end 1:end-1],:) )/3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = compute_curvabs(a)

u = a([2:end 1],:) - a([end 1:end-1],:);
s = sqrt(sum(u.^2,2));
s = [0;cumsum(s)];
s = s/s(end); 
s = s(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = compute_normal(a)

u = a([2:end 1],:) - a([end 1:end-1],:);
u = u ./ repmat( sqrt(sum(u.^2,2)), [1 2] );
u = [-u(:,2) u(:,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = tooth_profile(t, a)

d = 1/2-a;
e = double(t<d) + double(t>=d & t<a+d).*(d+a-t)/(a) + ...
        double(t>=1-a).*(t-(1-a))/(a); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = normal_extrude(a)

n = a-a(:,[2:end 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = genpolygon_regular(p,n)

theta = linspace(0,2*pi,p+1)';
theta(end) = [];
A = [cos(theta), sin(theta)];
a = genpolygon(A,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = genpolygon(A,n)

p = size(A,1);

a = [];
for i=1:p
    j = mod(i,p)+1;    
    if i<p
        k = floor(n/p);
    else
        k = n - (p-1)*floor(n/p);
    end
    t = (0:k-1)'/k;
    a(end+1:end+k,:) = (1-t)*A(i,:) + t*A(j,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cart2gear(a)

n = size(a,1);
[theta,r] = cart2pol(a(:,1), a(:,2));
theta = theta - theta(1);
theta(theta<0) = theta(theta<0) + 2*pi;
theta(end+1) = 2*pi; r(end+1) = r(1);

theta0 = (0:n-1)'/n*2*pi;
x = interp1(theta,r,theta0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = gear2cart(x)

n = length(x);
theta = (0:n-1)'/n*2*pi;
a = [x.*cos(theta), x.*sin(theta)];
