%%
% display in colors of complex polynomials

addpath('../toolbox/');
rep = MkResRep('animate');

if not(exist('cnt'))
    cnt = 1;
end

% grid
r = 1.5;
n = 501;
t = linspace(-r,r,n);
[Y,X] = meshgrid(t,t);
Z = X+1i*Y;


Proots = {};
for k=1:2
    Proots{k} = [];
    F = ones(n);
    clf; hold on;
    while true
        axis equal; axis([-r r -r r]);
        if length(Proots{k})>0
            PolyDisp(F,t);
        end
        if k==2
            plot(real(Proots{1}), imag(Proots{1}), 'b.', 'MarkerSize', 25);
        end
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if button==3
            break;
        end
        Proots{k}(end+1) = a+1i*b;
        F = F .* (Z-(a+1i*b));
    end
    Pcoefs{k} = poly(Proots{k});
end

% pad with zero
d = max(length(Pcoefs{1}), length(Pcoefs{2}));
for k=1:2
    Pcoefs{k} = [zeros(d-length(Pcoefs{k}),1); Pcoefs{k}(:)];
end
    
q = 50;
for i=1:q
    m = (i-1)/(q-1);
    P = (1-m)*Pcoefs{1} + m*Pcoefs{2};
    Pr = roots(P);
    F = polyval(P,Z);
    %
    clf;
    PolyDisp(F,t);
    plot(real(Pr), imag(Pr), 'b.', 'MarkerSize', 30);
    axis([-r r -r r]);
    drawnow;
    saveas(gcf, [rep 'anim-' num2str(cnt) '-' znum2str(i,2) '.png']);
end
cnt = cnt+1;

% create gifs
% AutoCrop(rep, ['anim-' num2str(cnt) '-'])
% Use % convert interp-*.png interp.gif to generate the gif using imagemagik


