function plot_neuron(x,y,w)

r = min(abs(w),1)*12;
if w>0
    plot(x, y, 'r', 'LineWidth', r);
elseif w<0
    plot(x, y, 'b', 'LineWidth', r);    
end
plot(x, y, 'k.', 'MarkerSize', 45);

end
