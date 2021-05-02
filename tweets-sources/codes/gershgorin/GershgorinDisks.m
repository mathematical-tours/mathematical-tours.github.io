%%
% Display Gershgorin circles
addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


n = 6;
h = 7;
for i=1:2
    A0{i} = ( randn(n) + 1i*randn(n) ) + h * diag( randn(n,1) + 1i*randn(n,1) );
end

M = max( [abs(eig(A0{1})); abs(eig(A0{2}))] );

c = exp(2i*pi*linspace(0,1,100));

q = 75;
for it=1:q
    t = (it-1)/(q-1);
    A = (1-t)*A0{1} + t*A0{2};
    clf; hold on;
    plot(eig(A), '.r', 'MarkerSize', 25);
    R = [];
    for i=1:n
        R(i) = sum( abs(A(i,[1:i-1,i+1:n])) );
        plot( A(i,i), 'b.', 'MarkerSize', 15);
        % fill( A(i,i) + R(i)*c, 'b', 'LineWidth', 2);
        fill( real(A(i,i) + R(i)*c), imag(A(i,i) + R(i)*c), 'b',  'EdgeColor', [.5 .5 1], 'LineWidth', 1, 'FaceAlpha', .1);
    end
    % plot the line
    for i=1:n
        plot( [A0{1}(i,i) A0{2}(i,i)], '-', 'Color', [1 1 1]*.5 );
    end
    
    scatter(real(diag(A)),imag(diag(A)), R, 'filled', 'b');
    
    axis equal;
    axis([-1 1 -1 1]*M*1.2);
    set(gca, 'XTick', [], 'YTick', []); box on;
    drawnow;
    mysaveas(it);
end