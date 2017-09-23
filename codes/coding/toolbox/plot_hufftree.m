function plot_hufftree(T,offs,S)

% plot_hufftree - plot a huffman tree 
%
%   plot_hufftree(T);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<2
    offs=0;
end
if nargin<3
    S = [];
end

hold on;
plot_tree(T{1},[0,0],1, offs, S);
hold off;
axis tight;
axis off;

end

%%
function plot_tree(T,x,j, offs, S)

tw = 15;
lw = 1.5;
ms = 20;

if not(iscell(T))
    str = num2str(T+offs);
    if not(isempty(S))
        str=S(T);
    end
    de = [-.02 -.2];
    % display a string
    u = text(x(1)+de(1),x(2)+de(2),str);
    set(u,'FontSize',tw);
else
    % display a split
    h0 = 1/2^j;
    h = 1 * h0;
    u = plot( [x(1) x(1)-h], [x(2) x(2)-1], 'k.-' );
    set(u,'MarkerSize',ms);
    set(u,'LineWidth',lw);
    u = plot( [x(1) x(1)+h], [x(2) x(2)-1], 'k.-' );
    set(u,'MarkerSize',ms);
    set(u,'LineWidth',lw);
    plot_tree(T{1},[x(1)-h x(2)-1],j+1, offs, S);
    plot_tree(T{2},[x(1)+h x(2)-1],j+1, offs, S);
end

end