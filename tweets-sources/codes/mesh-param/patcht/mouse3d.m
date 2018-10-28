function mouse3d(varargin)
% This function MOUSE3D enables mouse camera control on an certain figure
% axes. 
%
% Enable mouse control with mouse3d(axis-handle) or just mouse3d 
%
%
% MouseButtons
%  Left : Rotate
%  Right : Zoom
%  Center : Pan
% Keys
%  'r' : Change mouse rotation from inplane to outplane
%  'i' : Go back to initial view
%
% Example,
%   [X,Y,Z] = peaks(30);
%   surf(X,Y,Z)
%   colormap hsv
%   % Enable mouse control
%   mouse3d
%
% Function is written by D.Kroon University of Twente (July 2010)

if(nargin<1)
    handle=gca;
else
    handle=varargin{1};
    if(ishandle(handle))
        if(~strcmpi(get(handle,'Type'),'axes'))
            error('mouse3d:input','no valid axis handle');
        end
    else
        error('mouse3d:input','no valid axis handle');
    end
end

handles.figure1=get(handle,'Parent');
handles.axes1=handle;


mouse3d_OpeningFcn(gcf, handles, varargin);

% --- Executes just before mouse3d is made visible.
function mouse3d_OpeningFcn(hObject, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mouse3d (see VARARGIN)

% Choose default command line output for mouse3d
handles.output = hObject;

% UIWAIT makes mouse3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);
data.mouse_position_pressed=[0 0];
data.mouse_position=[0 0];
data.mouse_position_last=[0 0];
data.mouse_pressed=false;
data.mouse_button='';
data.firsttime=true;
data.mouse_rotate=true;
data=loadmousepointershapes(data);
data.handles=handles;
data.trans=[0 0 0];
data.Mview=[1 0 0 0;0 1 0 0; 0 0 1 0; 0 0 0 1];
setMyData(data);

set(data.handles.axes1,'ButtonDownFcn',@axes1_ButtonDownFcn);
set(data.handles.figure1,'WindowButtonMotionFcn',@figure1_WindowButtonMotionFcn);
set(data.handles.figure1,'WindowButtonUpFcn',@figure1_WindowButtonUpFcn);
set(data.handles.figure1,'KeyPressFcn',@figure1_KeyPressFcn);

function setViewMatrix()
data=getMyData; if(isempty(data)), return, end
Mview=data.Mview;
UpVector=Mview(1,1:3); 
Camtar=[0 0 0];
XYZ=Mview(2,1:3);
Forward=Mview(3,1:3);
    
UpVector=cross(cross(UpVector,Forward),Forward);
trans2=Mview(1:3,1:3)\data.trans(:);

Camtar=Camtar-Mview(1:3,4)'+trans2(1:3)'*data.scale;
XYZ=(XYZ-Mview(1:3,4)'+trans2(1:3)')*data.scale;
set(data.handles.axes1,'CameraUpVector', UpVector);
set(data.handles.axes1,'CameraPosition', XYZ+data.center);
set(data.handles.axes1,'CameraTarget', Camtar+data.center);
drawnow;


function setWindow()
data=getMyData; if(isempty(data)), return, end
c=get(data.handles.axes1,'Children');
jx=0; jy=0; jz=0;
if(~isempty(c))
    pmin=zeros(length(c),3);
    pmax=zeros(length(c),3);
    for i=1:length(c)
        c2=get(c(i));
        if(isfield(c2,'XData')),
            xd=c2.XData(:); xd(isnan(xd))=[];
            jx=jx+1;
            pmin(jx,1)= min(xd);
            pmax(jx,1)= max(xd);
        end
        if(isfield(c2,'YData')), 
            yd=c2.YData(:);  yd(isnan(yd))=[];
            jy=jy+1;
            pmin(jy,2)= min(yd);
            pmax(jy,2)= max(yd);
        end
        if(isfield(c2,'ZData')),
            zd=c2.ZData(:); zd(isnan(zd))=[];
            jz=jz+1;
            pmin(jz,3)= min(zd);
            pmax(jz,3)= max(zd);
        end
    end
    jx(jx==0)=1; jy(jy==0)=1; jz(jz==0)=1;
    pmin=[min(pmin(1:jx,1)) min(pmin(1:jy,2)) min(pmin(1:jz,3))];
    pmax=[max(pmax(1:jx,1)) max(pmax(1:jy,2)) max(pmax(1:jz,3))];
    data.center=(pmax+pmin)/2;
    data.scale=max(pmax-pmin)/2;
else
    data.center=[0 0 0];
    data.scale=1;
end
setMyData(data);

axis([-1 1 -1 1 -1 1]*data.scale+[data.center(1) data.center(1) data.center(2) data.center(2) data.center(3) data.center(3)]);
drawnow

set(data.handles.axes1,'CameraPositionMode','manual');
set(data.handles.axes1,'CameraUpVectorMode','manual');
set(data.handles.axes1,'CameraTargetMode','manual');
set(data.handles.axes1,'CameraViewAngleMode','manual');
set(data.handles.axes1,'PlotBoxAspectRatioMode','manual');
set(data.handles.axes1,'DataAspectRatioMode','manual');
set(data.handles.axes1,'CameraViewAngle',100);
set(get(data.handles.axes1,'Children'),'ButtonDownFcn',@axes1_ButtonDownFcn);
set(data.handles.axes1,'ButtonDownFcn',@axes1_ButtonDownFcn);
setViewMatrix()

function figure1_WindowButtonMotionFcn(hObject, eventdata)
cursor_position_in_axes();
data=getMyData(); if(isempty(data)), return, end
if(data.firsttime)
    data.firsttime=false; setMyData(data);
    setWindow();   
end

if(data.mouse_pressed)
    t1=(data.mouse_position_last(1)-data.mouse_position(1));
    t2=(data.mouse_position_last(2)-data.mouse_position(2));
    switch(data.mouse_button)
    case 'rotate1'
        R=RotationMatrix([t1 0 t2]);
        data.Mview=R*data.Mview;
        setMyData(data);
        setViewMatrix()
    case 'rotate2'
        R=RotationMatrix([0 0.5*(t1+t2) 0]);
        data.Mview=R*data.Mview;
        setMyData(data);
        setViewMatrix()
    case 'pan'
        data.trans=data.trans+[-t2/100 0 -t1/100];
        setMyData(data);
        setViewMatrix()
    case 'zoom'
        z=1-t2/100;
        R=ResizeMatrix([z z z]); 
        data.Mview=R*data.Mview;
        setMyData(data);
        setViewMatrix()
    otherwise
    end
end

function R=RotationMatrix(r)
% Determine the rotation matrix (View matrix) for rotation angles xyz ...
    Rx=[1 0 0 0; 0 cosd(r(1)) -sind(r(1)) 0; 0 sind(r(1)) cosd(r(1)) 0; 0 0 0 1];
    Ry=[cosd(r(2)) 0 sind(r(2)) 0; 0 1 0 0; -sind(r(2)) 0 cosd(r(2)) 0; 0 0 0 1];
    Rz=[cosd(r(3)) -sind(r(3)) 0 0; sind(r(3)) cosd(r(3)) 0 0; 0 0 1 0; 0 0 0 1];
    R=Rx*Ry*Rz;
    
function M=ResizeMatrix(s)
	M=[1/s(1) 0 0 0;
	   0 1/s(2) 0 0;
	   0 0 1/s(3) 0;
	   0 0 0 1];

function axes1_ButtonDownFcn(hObject, eventdata)
data=getMyData(); if(isempty(data)), return, end
handles=data.handles;
data.mouse_pressed=true;
data.mouse_button=get(handles.figure1,'SelectionType');
data.mouse_position_pressed=data.mouse_position;
if(strcmp(data.mouse_button,'normal'))
    if(data.mouse_rotate)
        data.mouse_button='rotate1';
        set_mouse_shape('rotate1',data);
    else
        data.mouse_button='rotate2';
        set_mouse_shape('rotate2',data);       
    end
end
if(strcmp(data.mouse_button,'open'))
end
if(strcmp(data.mouse_button,'extend'))
    data.mouse_button='pan';
    set_mouse_shape('pan',data);
end
if(strcmp(data.mouse_button,'alt'))
    data.mouse_button='zoom';
    set_mouse_shape('zoom',data);
end
setMyData(data);


function data=loadmousepointershapes(data)
I=[0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0;
 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 1 1 1 1 0 0 0 0 0 0;
 0 1 1 1 1 1 0 1 0 0 1 1 1 1 1 1; 1 1 1 1 0 0 0 1 0 0 0 0 0 0 1 1;
 1 1 1 1 0 0 0 1 0 0 0 0 0 1 1 1; 1 0 0 0 0 0 0 0 1 0 1 0 0 0 1 1;
 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 1; 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0];
I(I==0)=NaN; data.icon_mouse_rotate1=I;
I=[1 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0; 1 1 0 1 1 1 1 1 1 1 1 1 1 0 0 0;
 1 1 1 1 1 1 0 0 0 0 0 1 1 1 0 0; 1 0 0 1 1 0 0 0 0 0 0 0 0 1 0 0;
 1 0 0 1 1 0 0 0 0 0 0 0 0 1 0 0; 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0;
 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1;
 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 1; 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 1;
 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 1; 0 0 1 1 1 0 0 0 0 0 1 1 1 1 1 1;
 0 0 0 1 1 1 1 1 1 1 1 1 1 0 1 1; 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 1]; 
I(I==0)=NaN; data.icon_mouse_rotate2=I;
I=[0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0;
 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0; 1 0 0 1 1 1 1 1 0 0 1 0 0 0 0 0;
 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0; 1 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0 0 1 0 0 0 0 1 1 1 1 0 0 0 0 0;
 0 0 0 1 1 1 1 0 0 1 1 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1];
I(I==0)=NaN; data.icon_mouse_zoom=I;
I=[0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 1 1 0 1 1 0 0 0 0 0 0; 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 1 1 0 0 0 1 0 0 0 1 1 0 0 0;
 0 1 1 0 0 0 0 1 0 0 0 0 1 1 0 0; 1 0 0 1 1 1 1 1 1 1 1 1 0 0 1 0;
 0 1 1 0 0 0 0 1 0 0 0 0 1 1 0 0; 0 0 1 1 0 0 0 1 0 0 0 1 1 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 1 0 1 0 0 0 0 0 0;
 0 0 0 0 0 1 1 0 1 1 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
I(I==0)=NaN; data.icon_mouse_pan=I;

function set_mouse_shape(type,data)

switch(type)
case 'rotate1'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_rotate1,'PointerShapeHotSpot',round(size(data.icon_mouse_rotate1)/2))
    set(data.handles.figure1,'Pointer','custom');
case 'rotate2'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_rotate2,'PointerShapeHotSpot',round(size(data.icon_mouse_rotate2)/2))
    set(data.handles.figure1,'Pointer','custom');
case 'select_distance'
    set(data.handles.figure1,'Pointer','crosshair')
case 'select_landmark'
    set(data.handles.figure1,'Pointer','crosshair')
case 'select_roi'
    set(data.handles.figure1,'Pointer','crosshair')
case 'normal'
    set(data.handles.figure1,'Pointer','arrow')
case 'alt'
    set(data.handles.figure1,'Pointer','arrow')
case 'open'
    set(data.handles.figure1,'Pointer','arrow')
case 'zoom'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_zoom,'PointerShapeHotSpot',round(size(data.icon_mouse_zoom)/2))
    set(data.handles.figure1,'Pointer','custom');
case 'pan'
    set(gcf,'Pointer','custom','PointerShapeCData',data.icon_mouse_pan,'PointerShapeHotSpot',round(size(data.icon_mouse_pan)/2))
    set(data.handles.figure1,'Pointer','custom');
    
otherwise
    set(data.handles.figure1,'Pointer',type);
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata)
data=getMyData(); if(isempty(data)), return, end
if(data.mouse_pressed)
    data.mouse_pressed=false;
    setMyData(data);
end
set_mouse_shape('arrow',data)
    
function cursor_position_in_axes()
data=getMyData(); if(isempty(data)), return, end;
data.mouse_position_last=data.mouse_position;
p = get(0, 'PointerLocation');
data.mouse_position=[p(1, 1) p(1, 2)];
setMyData(data);

function setMyData(data)
% Store data struct in figure
setappdata(gcf,'data3d',data);

function data=getMyData()
% Get data struct stored in figure
data=getappdata(gcf,'data3d');

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata)
data=getappdata(gcf,'data3d');
switch(eventdata.Character)
    case 'i'
        data.trans=[0 0 0];
        data.Mview=[1 0 0 0;0 1 0 0; 0 0 1 0; 0 0 0 1];
        setMyData(data);
        setViewMatrix();
    case 'r'
        data.mouse_rotate=~data.mouse_rotate;
        setMyData(data)
    otherwise
end

