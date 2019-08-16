%%
% n-players matrix games.

addpath('../toolbox/');
rep = MkResRep();

n = 3; % #players
% A should be assymetric for 0 sum games
A = [0 1 1; -1 0 1; -1 -1 0];

% x_i -> min_x xi * A(i,:)*x
% F(x) = A*x + diag(A)*x

x0 = rand(3,1); 
x0 = [1 1 0]' * .4;

tau = .1;

niter = 50;
q = 50; % for anim
ndisp = round(linspace(1,niter,q));
kdisp = 1;
       
% display quantized colormap
Xgrad = x0;
Xextr = x0;


for it=1:niter
    % 'explicit'
    Xgrad(:,end+1) = Xgrad(:,end) - tau*A*Xgrad(:,end);
    % 'extra'
    x1 = Xextr(:,end) - tau*A*Xextr(:,end);
    Xextr(:,end+1) = Xextr(:,end) - tau*A*x1;

    if ndisp(kdisp)==it
        % display
        r = 15; % #levellines
        clf; hold on;
%         imagesc(t,t,F');
%         contour(t,t,F',linspace(min(F(:)), max(F(:)),r), 'k');
%         colormap(parula(r-1));
%         caxis([min(F(:)) max(F(:))]);
%         [x,y] = meshgrid(-1:.15:1,-1:.15:1);
%         quiver(x,y,-y,x,'k');
%         axis image; axis off;
        % plot curves
        plot3(Xgrad(1,:), Xgrad(2,:), Xgrad(3,:), '.-', 'LineWidth', 2, 'color', 'r', 'MarkerSize', 15);
        plot3(Xextr(1,:), Xextr(2,:), Xextr(3,:), '.-', 'LineWidth', 2, 'color', 'b', 'MarkerSize', 15);
        axis equal; 
        % axis tight;
         axis([-1 1 -1 1 -1 1]);
        view(3); box on;
        drawnow;
        % saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
        kdisp = kdisp+1;
    end
    
end
