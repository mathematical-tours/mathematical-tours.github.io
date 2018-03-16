function sphere_disp(V,F,f,g,alpha,beta)

% display elevation map on the sphere

if nargin<4
    rho = .3;
end
m = max(abs(g(:)));
if m>0
    g = g/m;
end
V1 = alpha*V + beta*repmat(g(:)', [3 1]) .*V;
% V1 = rho*repmat(g(:)', [3 1]) .*V;

opt.face_vertex_color = rescale(f(:));
plot_mesh(V1, F, opt);
colormap jet(256);

end