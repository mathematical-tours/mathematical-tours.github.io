
d = 6; % dimension

addpath('../toolbox/');
rep = MkResRep(num2str(d));

[C,I] = GenTesseract(d);

[C0,I0] = GenTesseract(3);

% random rotation
randn('state', 123);
Q = qr(randn(d));

eta = 5;
switch d
    case 4
        eta = 3;
    case 5
        eta = 2;
    case 6
        eta = 1;
end

S = randn(d);  S = eta*(S-S')/2;

%
q = 80;
for i=1:q
    t = (i-1)/q;
%    U = expm(t*logm(Q));
    U = expm(t*S);
    C1 = U*C;    
    % display the cube
    col = [1 0 0];
    clf; hold on;
    bb = 1.8;
    plot3(C0(1,:)*bb,C0(2,:)*bb,C0(3,:)*bb, 'w.', 'MarkerSize', 1e-3);
    PlotCube(C1,I, col);
    view(3); 
    axis equal;
    drawnow;
    % 
    saveas(gcf, [rep 'tesser-' znum2str(i,2) '.png']);
end


% AutoCrop(rep, ['tesser-']); 