% Section 8.4.1, Boyd & Vandenberghe "Convex Optimization"
% Almir Mutapcic - 10/05
% (a figure is generated)
%
% Given a finite set of points x_i in R^2, we find the minimum volume
% ellipsoid (described by matrix A and vector b) that covers all of
% the points by solving the optimization problem:
%
%           maximize     log det A
%           subject to   || A x_i + b || <= 1   for all i
%
% CVX cannot yet handle the logdet function, but this problem can be
% represented in an equivalent way as follows:
%
%           maximize     det(A)^(1/n)
%           subject to   || A x_i + b || <= 1   for all i
%
% The expression det(A)^(1/n) is SDP-representable, and is implemented
% by the MATLAB function det_rootn().


q = 40;

x = []; clf;
for it=1:q

    axis equal; axis([0 1 0 1]);
    box on; set(gca, 'Xtick', [], 'Ytick', []);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);
    if button==3
        break;
    end
    x(:,end+1) = [a;b];
    
    if it>2        
        % Create and solve the model
        cvx_begin quiet
        % cvx_precision best
        variable A(n,n) symmetric
        variable b(n)
        % maximize( det_rootn( A ) )
        maximize(log_det(A))
        subject to
        % norms( A * x + b * ones( 1, m ), 2 ) <= 1;
        % norms( A * x + b * ones( 1, m ) ) <= 1;
        for i=1:size(x,2)
            norm( A * x(:,i) + b,2 ) <= 1;
        end
        cvx_end
        % sum( (A*x+b).^2 )
        % Plot the results
        noangles = 200;
        angles   = linspace( 0, 2 * pi, noangles );
        z  = A \ [ cos(angles) - b(1) ; sin(angles) - b(2) ];
        clf; hold on;
        plot( x(1,:), x(2,:), 'r.', 'MarkerSize', 25 );
        plot( z(1,:), z(2,:), 'b-', 'LineWidth', 2 );
        axis equal; axis([0 1 0 1]);
        box on; set(gca, 'Xtick', [], 'Ytick', []);
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(it,q) '.png']);
    end
end