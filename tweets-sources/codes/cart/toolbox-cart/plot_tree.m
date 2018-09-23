function plot_tree(Tree)

% plot_tree - display a tree
%
%   plot_tree(Tree);
%
%   Copyright (c) 2007 Gabriel Peyre

J = length(Tree);
% branching factor
q = length(Tree{2})/length(Tree{1});

% edge
edgecolor = 'b';
% leaf
leafcolor = 'r.';
leafsize = 20;
% node
nodecolor = 'b.';
nodesize = 15;

delta = [-0.018,-0.045];
textsize = 20;

% clf;
hold on;
for j=1:J
    t = Tree{j};
    hj = 1-(j-1)*1/J;
    Hj = 1-j*1/J;
    nj = q^(j-1); % number of nodes
    xj = linspace(0,1,nj+2); xj(1) = []; xj(end) = [];
    Xj = linspace(0,1,nj*q+2); Xj(1) = []; Xj(end) = [];
    for i=1:nj        
        if t(i)==0
            % plot two branch
            for s=1:q
                plot([xj(i) Xj(q*(i-1)+s)], [hj Hj], edgecolor);
            end
%            plot([xj(i) Xj(2*i)], [hj   Hj], edgecolor);
            plot_point(xj(i), hj, nodecolor, nodesize);
            % plot the choice
            if 0
                h = text( xj(i)+delta(1), hj+delta(2), num2str(t(i)) );  
                set(h, 'FontSize', textsize);
                set(h, 'FontWeight', 'bold');
            end
            if j==J
                plot_point(Xj(2*i-1), Hj, leafcolor, leafsize);
                plot_point(Xj(2*i), Hj, leafcolor, leafsize);                
            end
        elseif t(i)==+1
            plot_point(xj(i), hj, leafcolor, leafsize);            
        elseif t(i)==-1
%            plot_point(xj(i), hj, 'k', leafsize/2);
        end
    end
end
% plot invisible bouding box
h = plot([0 1 1 0], [0 0 1 1]);
set(h, 'LineStyle', 'none');
hold off;
axis tight; axis off;

%%%
function plot_point(x,y,c,s)
h = plot(x, y, c);
set(h, 'MarkerSize', s);