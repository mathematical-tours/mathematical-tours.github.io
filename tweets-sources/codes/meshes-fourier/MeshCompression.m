addpath('../toolbox/');
rep = MkResRep();

name = 'duck';
name = 'moomoo';
[V,F] = read_off([name '.off']);
options.name = name; 
n = size(V,2);

E = [F([1 2],:) F([2 3],:) F([3 1],:)];
p = size(E,2);

% adjacency
W = sparse( E(1,:), E(2,:), ones(p,1) );
W = max(W,W');
% Laplacian
D = spdiags(sum(W)', 0, n,n);
L = D-W;

nb = 150;
opts.disp = 0;
[U,S] = eigs(L,nb,'SM',opts);
S = diag(S);



ilist = round(linspace(3,nb, 6));
tau=2.2; % saturation for display
clf;
for i=1:length(ilist)
    v = real(U(:,ilist(i)));
    v = clamp( v/std(v),-tau,tau );
    options.face_vertex_color = v;
    subplot(2,3,i);
    plot_mesh(V,F,options);
    shading interp; camlight; axis tight;
    colormap jet(256);
end

thresh = @(x,t)1/2+1/pi*atan((x-t)/.3);

q = 80;
mlist = 3 + (nb-1)*linspace(0,1,q).^3;
% mlist = unique(mlist); 
q = length(mlist);
for it=1:q
    r = mlist(it);
    w = diag(1-thresh(1:nb,r));
    V1 = (  (V * U) * w )*U';
    opt.face_vertex_color = ones(n,1);
    opt.name = name;
    clf;
    plot_mesh(V1,F, opt);
    shading interp; 
    shading faceted; 
    camlight; axis tight;
    colormap gray(256);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,3), '.png']);
end
   
