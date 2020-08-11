
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% # grids
n = 23;
n = 27;
% sampling on grid
p = 601;


x = linspace(-1,1,n);
y = linspace(-1,1,p);
[Y,X] = meshgrid(y,x); X = X'; Y = Y';



Disk2Half = @(z)(z+1i)./(1i*z+1);
Half2Disk = @(z)(z-1i)./(-1i*z+1);

Moeb  = @(a,b,c,d,z)(a*z+b)./(c*z+d);
MoebC = @(w1,w2,k,z)Moeb( w1-k*w2,(k-1)*w1*w2, 1-k, k*w1-w2,z );

mymode = 'half';
mymode = 'disk';
mymode = 'plane';


% for disk/half
Phi = [0 .2*pi];
B = [0 .9*exp(.3*2i*pi)];
% for plane
K = [1 .8*exp(2i*pi*0.6)];
W1 = [-.5,-.5];
W2 = [.5,.5];
% Grids
switch mymode
    case 'plane'
        G = X+1i*Y;
        H = Y+1i*X;
    case 'disk'
        G = exp(2i*pi*(X+1)/2).*(Y+1)/2;
        H = exp(2i*pi*(Y+1)/2).*(X+1)/2;
    case 'half'
        G = exp(2i*pi*(X+1)/2).*(Y+1)/2;
        H = exp(2i*pi*(Y+1)/2).*(X+1)/2;
        G = Disk2Half(G); 
        H = Disk2Half(H); 
end



q = 150;
lw=2; fs = 20;
damp = 0;
for it=1:q
    t = (it-1)/(q-1);
    switch mymode
        case 'plane'
            k = (1-t)*K(1) + t*K(2);
            w1 = (1-t)*W1(1) + t*W1(2);
            w2 = (1-t)*W2(1) + t*W2(2);
            f = @(z)MoebC(w1,w2,k,z)
        case 'disk'
            phi = (1-t)*Phi(1) + t*Phi(2);
            b = (1-t)*B(1) + t*B(2);
            m = exp(1i*phi);
            f = @(z)Moeb(m,m*b, conj(b), 1, z);  
        case 'half'
            phi = (1-t)*Phi(1) + t*Phi(2);
            b = (1-t)*B(1) + t*B(2);
            m = exp(1i*phi);
            f = @(z)Disk2Half( Moeb(m,m*b, conj(b), 1, Half2Disk(z)) );  
            
    end
    %
    clf; hold on;
    plot(f(G), 'color', [1 0 0]*(1-damp)+[1 1 1]*damp, 'LineWidth', lw);
    plot(f(H), 'color', [0 0 1]*(1-damp)+[1 1 1]*damp, 'LineWidth', lw);
    %    
    axis equal;  axis tight; 
    switch mymode
        case 'plane'
            plot(real([w1 w2]), imag([w1 w2]), 'k.', 'MarkerSize', 25);
            axis([-1 1 -1 1]*2);
        case 'disk'
            plot( exp(2i*pi*linspace(0,1,200)), 'k', 'LineWidth', 3 );
            axis([-1 1 -1 1]*1.2);
        case 'half'
            plot([-1 1]*2.5, [0 0], 'k', 'LineWidth', 3);
            axis([-2.5 2.5 -.2 5]);
    end    
    box on;
    set(gca, 'FontSize', fs); axis off;
    drawnow;
    mysaveas(it);
end
