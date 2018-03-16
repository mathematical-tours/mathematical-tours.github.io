rep = '../results/silouhette/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

name = 'elephant';
[V,F] = read_off([name '.off']);
m = size(F,2); n = size(V,2);

clf;
plot_mesh(V,F);
camlight;

% make a quantized rendeing

axis off;
saveas(gcf, 'tmp.png', 'png');
f = imread('tmp.png');

f = rescale(sum(f,3));
clf;
imagesc(f);
colormap gray(3);

% compute normal
[~,Nm] = compute_normal(V,F);
% compute center of faces
Fc = ( V(:,F(1,:)) + V(:,F(2,:)) + V(:,F(3,:)) )/3;
% inner product view->face / normals
u = get(gca, 'CameraPosition'); U = repmat(u(:),[1 m]);
IP = sum( (Fc-U).*Nm, 1 );

opt.face_vertex_color = double(IP>0)';
clf;
plot_mesh(V,F,opt);
camlight;


% compute edges wich are facing

e2f = compute_edge_face_ring(F);
[i,j,M] = find(e2f);
I = find(i<j);
i = i(I); j = j(I);
%
S = sign(IP(e2f(i,j))) .* sign(IP(e2f(j,i)));
I = find(S<0);
