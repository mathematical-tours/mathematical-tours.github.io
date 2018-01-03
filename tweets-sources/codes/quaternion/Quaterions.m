%%
% SLERP in quaterion spaces.


rep = '../results/quaternions/';
[~,~] = mkdir(rep);
if not(exist('test'))
    test = 0;
end
test = test+1;

q2rot = @(a,b,c,d)[...
        a^2+b^2-c^2-d^2, 2*b*c-2*a*d, 2*b*d+2*a*c; ...
        2*b*c+2*a*d, a^2-b^2+c^2-d^2, 2*c*d-2*a*b; ... 
        2*b*d-2*a*c, 2*c*d+2*a*b, a^2-b^2-c^2+d^2];
Q2Rot = @(q)q2rot(q(1),q(2),q(3),q(4));

normalize = @(z)z/norm(z);
q = normalize(randn(4,1));
U = Q2Rot(q); U*U';

Angle2Quat = @(u,t)[cos(t/2); u*sin(t/2)];

% slerp

dotp = @(u,v)sum(u.*v);
Slerp = @(u,v,t,h)(sin((1-t)*h)*u + sin(t*h) * v)/sin(h);
Slerp = @(u,v,t)Slerp(u,v,t,acos(dotp(u,v)));

u0 = normalize(randn(3,1)); t0 = -pi/2;
u1 = normalize(randn(3,1)); t1 = pi/3;
q0 = Angle2Quat(u0,t0);
q1 = Angle2Quat(u1,t1);

m = 16;
h = .8/m;
lw = 2;
clf; hold on;
for i=1:m
    t = (i-1)/(m-1);
    % interp
    qt = Slerp(q0,q1,t);
    Ut = Q2Rot(qt);
    % draw
    c = [t;0;0]; % center
    col = {'r' 'g' 'b'};
    for s=1:3
        d = c + Ut(:,s)*h;
        plot3([c(1) d(1)], [c(2) d(2)], [c(3) d(3)], col{s}, 'LineWidth', lw);
    end
end
axis equal;
axis off;
saveas(gcf, [rep 'quat-' num2str(test) '.eps'], 'epsc');



