
% create save repertory
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

% A Dubin's curve is a nearly kinemetically feasible path solution for car-like platform. The method explicitly find the trajectory composed of 3 segment: two curves and one staight line, or three curves. The curves are part of the circle. There are only 6 kind of the composition that proved to be minimun in length, and thease 6 types are called Dubin's curve.

r = .1;


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);
k = 4*2; x = rand(k,1)+1i*rand(k,1);
eta = .02; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);
theta = rand(k,1)*2*pi; thetav = .06 * (-1).^(1:k)';
q = 120;
col = distinguishable_colors(k+2);  col(4,:) = [];
for it=1:q
    % com
    clf; hold on;
    for j=1:k/2        
        p1 = [real(x(j)) imag(x(j)) theta(j)];
        p2 = [real(x(j+k/2)) imag(x(j+k/2)) theta(j+k/2)];
        q = dubins_curve(p1, p2, r, .01/5, true);
        q(end+1,:) = p2';
        % dusplay
        plot(q(:,1), q(:,2), 'color', col(j,:), 'LineWidth', 2);
        s = .1;
        plot([p1(1),p1(1) + s*cos(p1(3))],[p1(2),p1(2) + s*sin(p1(3))], 'color', col(j,:)*.7, 'LineWidth', 2);
        plot([p2(1),p2(1) - s*cos(p2(3))],[p2(2),p2(2) - s*sin(p2(3))], 'color', col(j,:)*.7, 'LineWidth', 2);
        plot(p1(1),p1(2), '.', 'color', col(j,:), 'MarkerSize', 20);
        plot(p2(1),p2(2), '.', 'color', col(j,:), 'MarkerSize', 20);
    end
    plot([0 0 1 1 0], [0 1 1 0 0], 'k');
    axis equal; axis off; axis([-.3 1.3 -.3 1.3]); 
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    % update
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
    theta = theta + thetav;
    % mysaveas(it);
end


