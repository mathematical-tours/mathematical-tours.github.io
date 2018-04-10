addpath('../toolbox/');
rep = MkResRep();


% grid
r = 2;
n = 501;
t = linspace(-r,r,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;

% click on points
z = [];
clf; hold on;
while true
    axis equal; axis([-r r -r r]);
    if length(z)>0
        DispOval(Z,z, t);
    end
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    z(end+1) = a+1i*b;
end