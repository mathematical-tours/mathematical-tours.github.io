function [P,Q,N,S] = LoadRoom(n)

% Walls of the cube
delta = .5/n; % small displacement to avoid singularity
t = linspace(delta,1-delta,n);
[B,A] = meshgrid(t,t); A = A(:); B = B(:); 
% small displ
t = linspace(0,1,n);
[B,A] = meshgrid(t,t); A = A(:); B = B(:); 
O = zeros(n*n,1); U = ones(n*n,1);

% the 6 walls, floor/left/right/top/bot/ceil, for display only
QX = [A; A; A; O; U; A];
QY = [B; O; U; A; A; B];
QZ = [O; B; B; B; B; U];
Q = cat(3, QX,QY,QZ);

% same but with small displacement
delta = .5/n; % small displacement to avoid singularity
A = .5 + (A-.5)*(1-1/n); B = .5 + (B-.5)*(1-1/n);
PX = [A; A; A; O; U; A];
PY = [B; O; U; A; A; B];
PZ = [O; B; B; B; B; U];
P = cat(3, PX,PY,PZ);

% normals to the walls
NX = [O; O; O; U;-U; O]; 
NY = [O; U;-U; O; O; O];
NZ = [U; O; O; O; O;-U];
N = cat(3, NX,NY,NZ);


% light source, a square on ceiling
W = zeros(n);
vx = [.2 .4]; vy = [.3 .7];
vx = round(vx*n); vy = round(vy*n);
W(vx(1):vx(2), vy(1):vy(2)) = 1; 
S = [O; O; O; O; O; W(:)]; 

end