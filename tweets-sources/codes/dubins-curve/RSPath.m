classdef RSPath
    properties (Constant)
        Types = [   
            'L', 'R', 'L', 'N', 'N' ;           %1
            'R', 'L', 'R', 'N', 'N' ;           %2
            'L', 'R', 'L', 'R', 'N' ;           %3
            'R', 'L', 'R', 'L', 'N' ;           %4
            'L', 'R', 'S', 'L', 'N' ;           %5
            'R', 'L', 'S', 'R', 'N' ;           %6
            'L', 'S', 'R', 'L', 'N' ;           %7
            'R', 'S', 'L', 'R', 'N' ;           %8
            'L', 'R', 'S', 'R', 'N' ;           %9
            'R', 'L', 'S', 'L', 'N' ;           %10
            'R', 'S', 'R', 'L', 'N' ;           %11
            'L', 'S', 'L', 'R', 'N' ;           %12
            'L', 'S', 'R', 'N', 'N' ;           %13
            'R', 'S', 'L', 'N', 'N' ;           %14
            'L', 'S', 'L', 'N', 'N' ;           %15
            'R', 'S', 'R', 'N', 'N' ;           %16
            'L', 'R', 'S', 'L', 'R' ;           %17
            'R', 'L', 'S', 'R', 'L'             %18
            ];
    end
    properties
type = repmat('N',[1,5]);
% duplicate array copies, ie ['N','N','N','N','N ']
t = 0; % The following 5 variables represent the path distance of the corresponding operation mode in type
        u = 0;
        v = 0;
        w = 0;
        x = 0;
        totalLength = 0;
    end
    methods
function obj = RSPath(type,t,u,v,w,x)% constructor
            obj.type = type;
            obj.t = t;
            obj.u = u;
            obj.v = v;
            obj.w = w;
            obj.x = x;
            obj.totalLength = sum(abs([t,u,v,w,x]));
        end
    end
end