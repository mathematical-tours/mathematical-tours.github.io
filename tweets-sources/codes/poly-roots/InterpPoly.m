%%
% Interpolate coefficients and vizualize roots -- new version.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


d = 30; % degree max

% click and play
for k=1:2
    clf; hold on; 
    if k>1
        plot(z0{1}, 'ko');
    end
    z0{k} = [];
    for i=1:d
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        if button==3
            d = i-1; 
            break;
        end
        plot(a,b, '.', 'MarkerSize', 15);
        z0{k}(end+1) = a+1i*b;
    end
    z0{k} = z0{k}(:);
    P0{k} = poly(z0{k});
end

Q = 500;
Z = [];
z = roots(P0{1});
for it=1:Q
    s = (it-1)/(Q-1);
    z1 = roots(P0{1}*(1-s) + P0{2}*s); z1 = z1(:);
    [~,I] = min( abs(z1-transpose(z)) ); 
    z = z1(I);    
    Z(:,it) = z(:);
end

Col = distinguishable_colors(d);
q = 70;
for it=1:q
    s = (it-1)/(q-1);
    i = round(1 + s*(Q-1));
    clf; hold on;
    for k=1:d
        plot(Z(k,:), '-', 'color', Col(k,:), 'LineWidth', 2);
        plot(Z(k,i), '.','color', Col(k,:), 'MarkerSize', 25); % 'color', [s 0 1-s], 
    end
    %
    axis equal; box on;
    e = .02;
    axis([-e 1+e -e 1+e]);
    set(gca, 'XTick', [], 'YTick', []); 
    drawnow;
    mysaveas(it);
end
