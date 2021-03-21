function PlotPath(path,veh)
    rmin = veh.MIN_CIRCLE;
    type = path.type;
    x = [];
    y = [];
    angle=[];
    seg = [path.t,path.u,path.v,path.w,path.x];
    pvec = [0,0,0];
    for i = 1:5        
        if type(i) == 'S'
            theta = pvec(3);
            dl = rmin*seg(i);
            dvec = [dl*cos(theta), dl*sin(theta), 0];
            dx = pvec(1)+linspace(0,dvec(1));
            dy = pvec(2)+linspace(0,dvec(2));
            x = [x,dx];
            y = [y,dy];
            pvec = pvec+dvec;
        elseif type(i) == 'L'
            theta = pvec(3);
            dtheta = seg(i);
            cenx = pvec(1)-rmin*sin(theta);
            ceny = pvec(2)+rmin*cos(theta);
            t = theta-pi/2+linspace(0,dtheta);
            dx = cenx+rmin*cos(t);
            dy = ceny+rmin*sin(t);
            x = [x,dx];
            y = [y,dy];
            angle=[angle,t];
            theta = theta+dtheta;
            pvec = [dx(end),dy(end),theta];
            dl = dtheta;
        elseif type(i) == 'R'
            theta = pvec(3);
            dtheta = -seg(i);
            cenx = pvec(1)+rmin*sin(theta);
            ceny = pvec(2)-rmin*cos(theta);
            t = theta+pi/2+linspace(0,dtheta);
            dx = cenx+rmin*cos(t);
            dy = ceny+rmin*sin(t);
            x = [x,dx];
            y = [y,dy];
            angle=[angle,t];
            theta = theta+dtheta;
            pvec = [dx(end),dy(end),theta];
            dl = -dtheta;
        else
            % do nothing
        end
        if dl > 0
            plot(dx,dy,'b');
        else
            plot(dx,dy,'r');
        end
        hold on
    end
    axis equal
    plot(0,0,'kx','LineWidth',2,'MarkerSize',10)
    plot(x(end),y(end),'ko', 'LineWidth',2,'MarkerSize',10)
%     veh = plot(x(1),y(1),'d','MarkerFaceColor','g','MarkerSize',10);
    videoFWriter = VideoWriter('Parking1.mp4','MPEG-4');
    open(videoFWriter);
[vehx,vehy] = getVehTran(x(1),y(1),angle(1));% Calculate the pose of the vehicle frame according to the pose of the center of the rear axle
h1 = plot(vehx,vehy,'r','LineWidth',4);% vehicle border
h2 = plot(x(1),y(1),'rx','MarkerSize',10);% vehicle rear axle center
    img = getframe(gcf);
    hold off
    pause(1)
    for k = 2:length(x)
        veh.XData = x(k);
        veh.YData = y(k);
        angle_2=angle(k);
        dl = norm([x(k)-x(k-1),y(k)-y(k-1)]);
        [vehx,vehy] = getVehTran(veh.XData,veh.YData,angle_2);
h1.XData = vehx;% Update the h1 image handle to add the x coordinates of the four corner points of the vehicle frame
        h1.YData = vehy;
h2.XData = veh.XData;% update h2 image handle, add the y coordinates of the four corner points of the vehicle frame
        h2.YData = veh.YData;
        writeVideo(videoFWriter,img);
        pause(dl)
    end    
end
