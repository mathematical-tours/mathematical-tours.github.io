%%
% Simple display of interpolation.

name = 'balanced';
name = 'unbalanced';
name = 'hellinger';

% create save repertory
addpath('../toolbox/');
rep = MkResRep(name);
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 1024; 
t = linspace(0,1,n);
s = .03; % width
G = @(m)exp(-(t-m).^2/(2*s^2));
GM = @(U,A)A(1)*G(U(1))+A(2)*G(U(2));
areac = @(f,c)area(t,f, 'FaceColor', c, 'EdgeColor', c, 'LineWidth', 1);


% first distrib
U = [.1 .75]; A = [1 .5];
V = [.25 .9]; B = [.7 .8];

q = 50; 
for it=1:q
    s = (it-1)/(q-1);
    
    % Unbalanced
    switch name
        case 'unbalanced'
            W = (1-s)*U + s*V;
            C = (1-s)*A + s*B;
            f = GM(W,C);
        case 'balanced'
            W = (1-s)*U + s*V;
            C = [.7 .5];
            f = GM(W,C) + .3 * G( (1-s)*.1 + s*.9 );
        case 'hellinger'
            f = (1-s)*GM(U,A) + s * GM(V,B);
    end
    
    clf; hold on;
    areac( GM(U,A), .8+.2*[1 0 0] );
    areac( GM(V,B), .8+.2*[0 0 1] );
    areac( f, [1-s 0 s] );
    axis([0 1 0 1.03]);
    set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 2/3 1]); box on;
    drawnow;
    mysaveas('anim',it); 

end

% AutoCrop(rep, 'anim-'); 