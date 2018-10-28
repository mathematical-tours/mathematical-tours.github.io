%%
% test for subdivision curves.

%% 
% Biblio:
% * [DynLevGre87]?N. Dyn, D. Levin and J.A. Gregory, <http://dx.doi.org/10.1016/0167-8396(87)90001-X _A 4-point interpolatory subdivision scheme for curve design_>, Computer Aided Geometric Design, 4(4), Pages 257-268, 1987.
% * [Chaikin74]?G. Chaikin, <http://dx.doi.org/10.1016/0146-664X(74)90028-8 _An algorithm for high speed curve generation_>. Computer Graphics and Image Processing, 3, 346-349, 1974.
% * [Riesen75]?R. Riesenfeld, <http://dx.doi.org/10.1016/0146-664X(75)90017-9 _On Chaikin's algorithm_>. Computer Graphics and Image Processing 4, 3, 304-310, 1975.
% * [DeslDub89]?G. Deslauriers and S. Dubuc. <http://dx.doi.org/10.1007/BF01889598 _Symmetric iterative interpolation processes_>. Constructive Approximation, 5(1):49-68, Dec. 1989.
% * [DynLevin02]?N. Dyn and D. Levin. <http://dx.doi.org/10.1017/S0962492902000028 _Subdivision schemes in geometric modelling_>. Acta Numerica, 11:73-144, Jan. 2002.
% * [DynGreLev91]?N. Dyn, J.A. Gregory, G. Levin, _Analysis of uniform binary subdivision schemes for curve design_, Constructive Approximation, 7(1), p. 127-147, 1991.


addpath('../toolbox/');

addpath('../toolbox/');
rep = MkResRep();


ms = 20; lw = 1.5;
myplot = @(f,c)plot(f([1:end 1]), c, 'LineWidth', lw, 'MarkerSize', ms);
myaxis = @(rho)axis([-rho 1+rho -rho 1+rho], 'off');




%%
% First we create a 2D closes polygon, which will be a cage used to perform
% 2D shape deformation.

clf; hold on;
f0 = [];
while true
    axis equal; axis([0 1 0 1]);  % axis off;
    [a,b,button] = ginput(1);  
    if button==3
        break;
    end
    plot(a,b, '.', 'MarkerSize', 25); 
    f0(end+1) = a + 1i*b;
    plot(f0, 'r', 'LineWidth', 2);
end
k = size(f0,2);
plot(f0([1:end,1]), 'r', 'LineWidth', 2);
f1 = [];
for it=1:k
    axis([0 1 0 1]); axis equal; axis off;
    [a,b,button] = ginput(1);
    plot(a,b, '*', 'MarkerSize', 25);   
    f1(end+1) = a + 1i * b;
end
f0 = f0(:); f1 = f1(:);


% basic step
subdivide = @(f,h)cconvol( upsampling(f), h);

% bi-cubic
name = 'quadratic';
name = 'cubic';
name = 'interpolating';
switch name
    case 'cubic'
        h = [1 4 6 4 1]; % cubic B-spline
        h = 2*h/sum(h(:));
    case 'quadratic'
        % [Chaikin74]
        % w=3 == quadratic B-spline
        hcc = @(w)[1 w w 1]/(1+w);
        h=hcc(3);
    case 'interpolating'
        % four-point interpolation kernel // DynLevGre87
        % w=1/16 == cubic B-spline interpolation.
        h4pt = @(w)[-w, 0, 1/2+w, 1, 1/2+w, 0, -w];
        h = h4pt(1/16);
end

q = 50;  % animation
jlist = [1 2 5]; 
f = f0;
for j=jlist
    for it=1:q
        t = (it-1)/(q-1);
        ft = f0*(1-t) + f1*t;
        
        f = ft;
        for isub=1:j
            f = subdivide(f,h);
        end
        
        st = '.-'; 
        if j==jlist(end)
           st = '-'; 
        end
        clf; hold on;
        plot(ft([1:end 1]), '.-', 'Color', [1 1 1]*.5, 'LineWidth', 1, 'MarkerSize', 25);
        plot(f([1:end 1]), st, 'Color', [t 0 1-t], 'LineWidth', 2, 'MarkerSize', 25);
        axis equal; myaxis(0.02);
        drawnow;
        saveas(gcf, [rep name '-' num2str(j) '-' znum2str(it,2) '.png'], 'png');
    end
end

AutoCrop(rep, name);

