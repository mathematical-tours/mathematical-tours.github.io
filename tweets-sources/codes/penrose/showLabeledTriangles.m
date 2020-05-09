function showLabeledTriangles(T)
%showLabeledTriangles Show triangles with type and side labels.
%
%   showLabeledTriangles(T) shows the outline of each triangle contained
%   in the input table. Each row of the input table has the form
%   returned by aTriangle, apTriangle, bTriangle, or bpTriangle.  Each
%   displayed triangled is labeled according to type: A, A', B, or B'.
%   Each triangle side is marked with a symbol that helps indicate a
%   correct or incorrect tiling.
%
%   EXAMPLE
%
%   Show how a B triangle is decomposed into three triangles (of type A,
%   B, and B') according to Penrose tiling rules.
%
%       t = bTriangle([],-1,1);
%       t1 = decomposeTriangles(t);
%       showLabeledTriangles(t1)

%   Copyright 2018 The MathWorks, Inc.

showTriangles(T);

vertices = [T.Left T.Apex T.Right];
centroids = mean(vertices,2);
cx = real(centroids);
cy = imag(centroids);

labels = cellstr(T.Type);
labels = strrep(labels,'Ap','A''');
labels = strrep(labels,'Bp','B''');

hold on
text(cx,cy,labels,...
   'HorizontalAlignment','center',...
   'VerticalAlignment','middle',...
   'FontWeight','bold',...
   'FontSize',14);

labelSides(T(T.Type == "A",:),...
    0.5,{'LineStyle','none','Marker','o','Color','k','MarkerSize',15},...
    0.5,{'LineStyle','none','Marker','s','Color','k','MarkerSize',12},...
    0.6,{'LineStyle','none','Marker','*','Color','k','MarkerSize',14});

labelSides(T(T.Type == "Ap",:),...
    0.5,{'LineStyle','none','Marker','s','Color','k','MarkerSize',12},...
    0.5,{'LineStyle','none','Marker','o','Color','k','MarkerSize',15},...
    0.4,{'LineStyle','none','Marker','*','Color','k','MarkerSize',14});
    
labelSides(T(T.Type == "B",:),...
    0.5,{'LineStyle','none','Marker','s','Color','k','MarkerSize',12},...
    0.5,{'LineStyle','none','Marker','o','Color','k','MarkerSize',15},...
    0.6,{'LineStyle','none','Marker','p','Color','k','MarkerSize',16});

labelSides(T(T.Type == "Bp",:),...
    0.5,{'LineStyle','none','Marker','o','Color','k','MarkerSize',15},...
    0.5,{'LineStyle','none','Marker','s','Color','k','MarkerSize',12},...
    0.4,{'LineStyle','none','Marker','p','Color','k','MarkerSize',16});

hold off

function labelSides(T,alpha_left,line_params_left,alpha_right,line_params_right,alpha_base,line_params_base)
side1_label_point = T.Apex + alpha_left*(T.Left - T.Apex);
plot(real(side1_label_point),imag(side1_label_point),line_params_left{:});

side2_label_point = T.Apex + alpha_right*(T.Right - T.Apex);
plot(real(side2_label_point),imag(side2_label_point),line_params_right{:});

base_label_point = T.Left + alpha_base*(T.Right - T.Left);
plot(real(base_label_point),imag(base_label_point),line_params_base{:});

