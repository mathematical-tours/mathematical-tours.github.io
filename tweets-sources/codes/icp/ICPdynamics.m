%% 
% display the dynamics of icp

% load two shapes
addpath('../toolbox/');
addpath('./mexEMD/');
rep = MkResRep();

match_mode = 'nn';
match_mode = 'oracle';
match_mode = 'ot';

name = 'cat';
name = 'bunny';
name = 'horse';
n0 = 500;
p = 700; % #points
f = load_image(name, n0);
f = rescale(sum(f,3)); f = f>.5;
if f(1)==1
    f = 1-f;
end
% random points
t = linspace(-1/2,1/2,n0);
[Y,X] = meshgrid(t,t);
I = find(f==1); I = I(randperm(length(I)));
I = I(1:p);
%
x = Y(I);
y = -X(I);

clf;
plot(x,y, '.', 'MarkerSize', 20);
axis equal; axis off;


% random rotation
[U0,~] = qr(randn(2));
U0 = U0 * diag([1 det(U0)]);
U0 = real( expm(.5*logm(U0)) ); % makes it smaller
v0 = .4*[1.7 1.3]';

Z = U0*[x,y]' + v0; 
x1 = Z(1,:)'; y1 = Z(2,:)';

tau = .45;
tau = .15;
q = 50;

xa = x; ya = y;
for it=1:q        
    % do registration
    D = abs( (xa+1i*ya) - transpose(x1+1i*y1) );
    switch match_mode
        case 'oracle'
            I = 1:p;
        case 'nn'
            [~,I] = min(D, [], 2);
        case 'ot'
            [cost,gamma] = mexEMD(ones(p,1)/p,ones(p,1)/p,D.^2);
            [J,K,gammaij] = find(gamma);
            Ji = zeros(p,1); Ji(J) = 1:p;
            I = K(Ji);
    end
    %
    s = (it-1)/(q-1);
    ms = 15;
    clf; hold on;
    g = .5;
    plot(x,y, '.', 'color', g*[1 1 1] + (1-g)*[0 0 1], 'MarkerSize', ms);
    plot(x1,y1, '.','color', g*[1 1 1] + (1-g)*[1 0 0], 'MarkerSize', ms);
    plot(xa,ya, '.', 'color', [s 0 1-s], 'MarkerSize', ms);
    % display a subset of OT connection
    S = round(linspace(1,p,100));
    for k=S
        plot([xa(k) x1(I(k))], [ya(k) y1(I(k))], '-', 'color', [1 1 1]*0);
    end        
    %
    axis equal; axis off;
    axis([-.5 1.2 -.5 .9]);
    axis([-.6 1.1 -.5 1.1]);
    set(gca, 'PlotBoxAspectRatio', [1 1 1]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
    
    % move points
    [U,v,x2,y2] = FitRotation(xa,ya,x1(I),y1(I));
    Ut = real( expm(tau*logm(U)) );    
    Z = Ut*[xa ya]' + tau*v;
    xa = Z(1,:)'; ya = Z(2,:)';
end
