%%
% Solve the radiosity equation inside a cube

addpath('../toolbox/');
rep = MkResRep();

n = 30; % discretization 
p = 6*n*n; % total #points

% helpers
dotp = @(u,v)sum(u.*v,3);

% load a square box room
[P,Q,N,S] = LoadRoom(n);

% sets of normalized ray vectors between pairs of points
pairw = @(a)repmat(a,[1 p])-repmat(permute(a, [2 1 3]),[p 1]);
R = pairw(P); 
D = max(1e-9, sqrt(sum(R.^2,3)));
R = R ./ repmat(D, [1 1 3] );

% conservation law kernel
N1 = repmat(N, [1 p]); N2 = repmat(permute(N, [2 1 3]), [p 1]);
K = max( -dotp(N1, R), 0) .* max( dotp(N2, R), 0) ./ (D).^2;
% ensure conservation of light
K = K ./ repmat( sum(K,2), [1 p] );

% reflectance material (should be <1 to ensure contractance)
rho = .9;

% positions
oc1 = { [.3 .2 .25] [.7 .2 .85]  }; 
oc2 = { [.7 .6 .8]  [.2 .8 .2] }; 
r = [.1 .15]; % radius

q = 40; % #frames
elev = 55;
rot = [40 130];
% delta for each frame in term of rotation
delta_rot = .5;

for i=1:q
    t = (i-1)/(q-1);
    oc = { (1-t)*oc1{1} + t*oc2{1}, (1-t)*oc1{2} + t*oc2{2} }; 
    % solve
    V = ComputeVisibilty(oc,r, P, R);
    K1 = K .* V;  K1 = K1 ./ repmat( sum(K1,2), [1 p] );
    f = (eye(p) - rho * K1 ) \ S; 
    % render
    clf; RenderScene(Q,f,S,oc,r);
    for j=1:length(rot)
        view(rot(j)+delta_rot*(i-1),elev);
        drawnow;
        saveas(gcf, [rep 'radiosity-vp' num2str(j) '-' znum2str(i,2) '.png' ]);
    end
end




if 0
% Taylor serie    
% f = (eye(p) - rho * K1 ) \ S; 
K = 3; 
f = zeros(p,1); 
for k=1:K
    f = S + rho*K1*f;
end
clf; RenderScene(Q,f,S,oc,r);

k = 0;
while true
    k = k+1;
    view(k,elev); axis([0 1 0 1 0 1]); drawnow;
end
end


% AutoCrop(rep, ['radiosity-'])
