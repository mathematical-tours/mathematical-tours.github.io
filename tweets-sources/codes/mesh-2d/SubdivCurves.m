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

rep = '../results/subdivision-curves/';
[~,~] = mkdir(rep);

ms = 20; lw = 1.5;
myplot = @(f,c)plot(f([1:end 1]), c, 'LineWidth', lw, 'MarkerSize', ms);
myaxis = @(rho)axis([-rho 1+rho -rho 1+rho], 'off');

% input coarse
f0 =    [0.11 0.18 0.26 0.36 0.59 0.64 0.80 0.89 0.58 0.22 0.18 0.30 0.58 0.43 0.42]' + ...
   1i * [0.91 0.55 0.91 0.58 0.78 0.51 0.81 0.56 0.10 0.16 0.35 0.42 0.40 0.24 0.31]';
f0 = rescale(real(f0),.01,.99) + 1i * rescale(imag(f0),.01,.99);

% basic step
subdivide = @(f,h)cconvol( upsampling(f), h);

% bi-cubic
name = 'interpolating';
name = 'quadratic';
name = 'cubic';
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


Jmax = 4;
f = f0;
for j=0:Jmax
    clf; hold on;
	plot(f0([1:end 1]), 'r.--', 'LineWidth', 1.5, 'MarkerSize', 25);
    plot(f([1:end 1]), 'k.-', 'LineWidth', 2, 'MarkerSize', 20);
    myaxis(0.02);
    saveas(gcf, [rep name '-' num2str(j) '.eps'], 'epsc');    
    f = subdivide(f,h);
end

