function plot_colored(B,alph,lw)

if nargin<2
    alph = .1;
end
if nargin<2
    lw = 1;
end
Q = size(B,1);
c1 = [1 0 0]; c0 = [0 0 1];
hold on;
for k=1:Q-1
    t = (k-1)/(Q-2);
    col = (1-t)*c0 + t*c1;
    u = real(B(k:k+1,:));
    v = imag(B(k:k+1,:));
%    plot( u,v, 'color', col, 'LineWidth', lw );    
	patch(u,v, 'r', 'EdgeColor',col, 'EdgeAlpha',alph,'FaceColor','none','Linesmoothing','on', 'LineWidth', lw);
end


end