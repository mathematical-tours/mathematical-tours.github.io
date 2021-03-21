function path = FindRSPath(x,y,phi,veh)
    rmin = veh.MIN_CIRCLE; %minimum turning radius
    x = x/rmin;
    y = y/rmin;
% traverse 5 methods to reach the target point, and then select the shortest path
    [isok1,path1] = CSC(x,y,phi);
    [isok2,path2] = CCC(x,y,phi);
    [isok3,path3] = CCCC(x,y,phi);
    [isok4,path4] = CCSC(x,y,phi);
    [isok5,path5] = CCSCC(x,y,phi);
    isoks = [isok1, isok2, isok3, isok4, isok5];
    paths = {path1, path2, path3, path4, path5};
    Lmin = inf;
% Find the curve with the shortest 5 paths
    for i = 1:5
        if isoks(i) == true
            elem = paths{i};
            if Lmin > elem.totalLength
                Lmin = elem.totalLength;
                path = elem;
            end
        end
    end
%     PlotPath(path,veh);
end

% The range of control angle x is [-pi,pi]
function v = mod2pi(x)
    v = rem(x,2*pi);
    if v < -pi
        v = v+2*pi;
    elseif v > pi
        v = v-2*pi;
    end
end

% formula 8.6
function [tau,omega] = tauOmega(u,v,xi,eta,phi)
    delta = mod2pi(u-v);
    A = sin(u)-sin(delta);
    B = cos(u)-cos(delta)-1;
    t1 = atan2(eta*A-xi*B,xi*A+eta*B);
    t2 = 2*(cos(delta)-2*cos(v)-2*cos(u))+3;
    if t2 < 0
        tau = mod2pi(t1+pi);
    else
        tau = mod2pi(t1);
    end
    omega = mod2pi(tau-u+v-phi);
end

% formula 8.1
function [isok,t,u,v] = LpSpLp(x,y,phi)
[t,u] = cart2pol(x-sin(phi),y-1+cos(phi));% Convert Cartesian coordinates to polar coordinates, return theta and rho, The paper returns [u, t], because the order of the values ??returned by the cart2pol function is different from the original text, the meaning of the variable is still the same, t stands for radians, u stands for the straight distance
if t >= 0% must turn left, t>=0 means turn left
        v = mod2pi(phi-t);
if v >= 0% symbol means forward and backward
            isok = true;
            return
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

% formula 8.2
function [isok,t,u,v] = LpSpRp(x,y,phi)
    [t1,u1] = cart2pol(x+sin(phi),y-1-cos(phi));
    if u1^2 >= 4
        u = sqrt(u1^2-4);
        theta = atan2(2,u);
        t = mod2pi(t1+theta);
        v = mod2pi(t-phi);
        if t >= 0 && v >= 0% symbol means forward and backward
            isok = true;
            return
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

function [isok,path] = CSC(x,y,phi)
    Lmin = inf;
    type = repmat('N',[1,5]);
    path = RSPath(type,0,0,0,0,0);
    [isok,t,u,v] = LpSpLp(x,y,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(15,:),t,u,v,0,0);
        end
    end
    [isok,t,u,v] = LpSpLp(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(15,:),-t,-u,-v,0,0);
        end
    end
    [isok,t,u,v] = LpSpLp(x,-y,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(16,:),t,u,v,0,0);
        end
    end
    [isok,t,u,v] = LpSpLp(-x,-y,phi); % timeflp + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(16,:),-t,-u,-v,0,0);
        end
    end
    [isok,t,u,v] = LpSpRp(x,y,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(13,:),t,u,v,0,0);
        end
    end
    [isok,t,u,v] = LpSpRp(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(13,:),-t,-u,-v,0,0);
        end
    end
    [isok,t,u,v] = LpSpRp(x,-y,-phi); % reflect 
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(14,:),t,u,v,0,0);
        end
    end
    [isok,t,u,v] = LpSpRp(-x,-y,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(14,:),-t,-u,-v,0,0);
        end
    end
    if Lmin == inf
        isok = false;
    else
        isok = true;
    end
end

% formula 8.3/8.4
function [isok,t,u,v] = LpRmL(x,y,phi)
    xi = x-sin(phi);
    eta = y-1+cos(phi);
    [theta,u1] = cart2pol(xi,eta);
    if u1 <= 4
        u = -2*asin(u1/4);
        t = mod2pi(theta+u/2+pi);
        v = mod2pi(phi-t+u);
        if t >= 0 && u <= 0
            isok = true;
            return
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

function [isok,path] = CCC(x,y,phi)
    Lmin = inf;
    type = repmat('N',[1,5]);
    path = RSPath(type,0,0,0,0,0);
    [isok,t,u,v] = LpRmL(x,y,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(1,:),t,u,v,0,0);
        end
    end
    [isok,t,u,v] = LpRmL(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(1,:),-t,-u,-v,0,0);
        end
    end
    [isok,t,u,v] = LpRmL(x,-y,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(2,:),t,u,v,0,0);
        end
    end
    [isok,t,u,v] = LpRmL(-x,-y,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(2,:),-t,-u,-v,0,0);
        end
    end
    % backwards
    xb = x*cos(phi)+y*sin(phi);
    yb = x*sin(phi)-y*cos(phi);
    [isok,t,u,v] = LpRmL(xb,yb,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(1,:),v,u,t,0,0);
        end
    end
    [isok,t,u,v] = LpRmL(-xb,yb,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(1,:),-v,-u,-t,0,0);
        end
    end
    [isok,t,u,v] = LpRmL(xb,-yb,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(2,:),v,u,t,0,0);
        end
    end
    [isok,t,u,v] = LpRmL(-xb,-yb,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(2,:),-v,-u,-t,0,0);
        end
    end
    if Lmin == inf
        isok = false;
    else
        isok = true;
    end
end

% formula 8.7,tauOmega() is formula 8.6
function [isok,t,u,v] = LpRupLumRm(x,y,phi)
    xi = x+sin(phi);
    eta = y-1-cos(phi);
    rho = (2+sqrt(xi^2+eta^2))/4;
    if rho <= 1
        u = acos(rho);
        [t,v] = tauOmega(u,-u,xi,eta,phi);
    if t >= 0 && v <= 0% The symbols represent forward and backward
            isok = true;
            return
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

% formula 8.8
function [isok,t,u,v] = LpRumLumRp(x,y,phi)
    xi = x+sin(phi);
    eta = y-1-cos(phi);
    rho = (20-xi^2-eta^2)/16;
    if rho >= 0 && rho <= 1
        u = -acos(rho);
        if u >= pi/2
            [t,v] = tauOmega(u,u,xi,eta,phi);
            if t >=0 && v >=0
                isok = true;
                return
            end
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

function [isok,path] = CCCC(x,y,phi)
    Lmin = inf;
    type = repmat('N',[1,5]);
    path = RSPath(type,0,0,0,0,0);
    [isok,t,u,v] = LpRupLumRm(x,y,phi);
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(3,:),t,u,-u,v,0);
        end
    end
    [isok,t,u,v] = LpRupLumRm(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(3,:),-t,-u,u,-v,0);
        end
    end
    [isok,t,u,v] = LpRupLumRm(x,-y,-phi); % reflect
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(4,:),t,u,-u,v,0);
        end
    end
    [isok,t,u,v] = LpRupLumRm(-x,-y,phi); % timeflip + reflect
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(4,:),-t,-u,u,-v,0);
        end
    end
    [isok,t,u,v] = LpRumLumRp(x,y,phi);
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(3,:),t,u,u,v,0);
        end
    end
    [isok,t,u,v] = LpRumLumRp(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(3,:),-t,-u,-u,-v,0);
        end
    end
    [isok,t,u,v] = LpRumLumRp(x,-y,-phi); % reflect
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(4,:),t,u,u,v,0);
        end
    end
    [isok,t,u,v] = LpRumLumRp(-x,-y,phi); % timeflip + reflect
    if isok
        L = abs(t)+2*abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(4,:),-t,-u,-u,-v,0);
        end
    end
    if Lmin == inf
        isok = false;
    else
        isok = true;
    end
end

% formula 8.9
function [isok,t,u,v] = LpRmSmLm(x,y,phi)
    xi = x-sin(phi);
    eta = y-1+cos(phi);
    [theta,rho] = cart2pol(xi,eta);
    if rho >= 2
        r = sqrt(rho^2-4);
        u = 2-r;
        t = mod2pi(theta+atan2(r,-2));
        v = mod2pi(phi-pi/2-t);
        if t >= 0 && u <= 0 && v <= 0
            isok = true;
            return
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

% formula 8.10
function [isok,t,u,v] = LpRmSmRm(x,y,phi)
    xi = x+sin(phi);
    eta = y-1-cos(phi);
    [theta,rho] = cart2pol(-eta,xi);
    if rho >= 2
        t = theta;
        u = 2-rho;
        v = mod2pi(t+pi/2-phi);
        if t >= 0 && u <= 0 && v <= 0
            isok = true;
            return
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

function [isok,path] = CCSC(x,y,phi)
    Lmin = inf;
    type = repmat('N',[1,5]);
    path = RSPath(type,0,0,0,0,0);
    [isok,t,u,v] = LpRmSmLm(x,y,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(5,:),t,-pi/2,u,v,0);
        end
    end
    [isok,t,u,v] = LpRmSmLm(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(5,:),-t,pi/2,-u,-v,0);
        end
    end
    [isok,t,u,v] = LpRmSmLm(x,-y,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(6,:),t,-pi/2,u,v,0);
        end
    end
    [isok,t,u,v] = LpRmSmLm(-x,-y,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(6,:),-t,pi/2,-u,-v,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(x,y,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(9,:),t,-pi/2,u,v,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(-x,y,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(9,:),-t,pi/2,-u,-v,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(x,-y,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(10,:),t,-pi/2,u,v,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(-x,-y,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(10,:),-t,pi/2,-u,-v,0);
        end
    end
    % backwards
    xb = x*cos(phi)+y*sin(phi);
    yb = x*sin(phi)-y*cos(phi);
    [isok,t,u,v] = LpRmSmLm(xb,yb,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(7,:),v,u,-pi/2,t,0);
        end
    end
    [isok,t,u,v] = LpRmSmLm(-xb,yb,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(7,:),-v,-u,pi/2,-t,0);
        end
    end
    [isok,t,u,v] = LpRmSmLm(xb,-yb,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(8,:),v,u,-pi/2,t,0);
        end
    end
    [isok,t,u,v] = LpRmSmLm(-xb,-yb,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(8,:),-v,-u,pi/2,-t,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(xb,yb,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(11,:),v,u,-pi/2,t,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(-xb,yb,-phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(11,:),-v,-u,pi/2,-t,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(xb,-yb,-phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(12,:),v,u,-pi/2,t,0);
        end
    end
    [isok,t,u,v] = LpRmSmRm(-xb,-yb,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(12,:),-v,-u,pi/2,-t,0);
        end
    end
    if Lmin == inf
        isok = false;
    else
        isok = true;
    end
end

% formula 8.11
function [isok,t,u,v] = LpRmSLmRp(x,y,phi)
    xi = x+sin(phi);
    eta = y-1-cos(phi);
    [~,rho] = cart2pol(xi,eta);
    if rho >= 2
        u = 4-sqrt(rho^2-4);
        if u <= 0
            t = mod2pi(atan2((4-u)*xi-2*eta,-2*xi+(u-4)*eta));
            v = mod2pi(t-phi);
            if t >= 0 && v >= 0
                isok = true;
                return
            end
        end
    end
    isok = false;
    t = 0;
    u = 0;
    v = 0;
end

function [isok,path] = CCSCC(x,y,phi)
    Lmin = inf;
    type = repmat('N',[1,5]);
    path = RSPath(type,0,0,0,0,0);
    [isok,t,u,v] = LpRmSLmRp(x,y,phi);
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(17,:),t,-pi/2,u,-pi/2,v);
        end
    end
    [isok,t,u,v] = LpRmSLmRp(x,y,phi); % timeflip
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(17,:),-t,pi/2,-u,pi/2,-v);
        end
    end
    [isok,t,u,v] = LpRmSLmRp(x,y,phi); % reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(18,:),t,-pi/2,u,-pi/2,v);
        end
    end
    [isok,t,u,v] = LpRmSLmRp(x,y,phi); % timeflip + reflect
    if isok
        L = abs(t)+abs(u)+abs(v);
        if Lmin > L
            Lmin = L;
            path = RSPath(RSPath.Types(18,:),-t,pi/2,-u,pi/2,-v);
        end
    end
    if Lmin == inf
        isok = false;
    else
        isok = true;
    end
end



% calculate the pose of the vehicle frame based on the pose of the center of the rear axle
function [x,y] = getVehTran(x,y,theta)
    W = 0.4;
    LF = 0.1;
    LB = 0.1;

% The frame of the vehicle is determined by four corner points
Cornerfl = [LF, W/2];% front left corner
Cornerfr = [LF, -W/2];% right front corner
Cornerrl = [-LB, W/2];% rear left corner
Cornerrr = [-LB, -W/2];% right rear corner point
Pos = [x,y];% center coordinates of rear axle
dcm = angle2dcm(-theta, 0, 0);% Calculate the rotation matrix of the four corner points, because it is part of a rigid body, the rotation matrix is ??the same, convert the angle to the direction cosine matrix, The rotation order is ZYX
    
tvec = dcm*[Cornerfl';0];% rotation transformation, column vector formed by Cornerfl rotation, position vector 3*1, the last one is z coordinate
    tvec = tvec';
Cornerfl = tvec(1:2)+Pos;% translation transformation
    
    tvec = dcm*[Cornerfr';0];
    tvec = tvec';
    Cornerfr = tvec(1:2)+Pos;
    
    tvec = dcm*[Cornerrl';0];
    tvec = tvec';
    Cornerrl = tvec(1:2)+Pos;
    
    tvec = dcm*[Cornerrr';0];
    tvec = tvec';
    Cornerrr = tvec(1:2)+Pos;
    
     
% returns the x,y coordinates of the four corner points of the vehicle frame
    x = [Cornerfl(1),Cornerfr(1),Cornerrr(1),Cornerrl(1),Cornerfl(1)];
    y = [Cornerfl(2),Cornerfr(2),Cornerrr(2),Cornerrl(2),Cornerfl(2)];
%     patch(x,y,[1 1 0])
end