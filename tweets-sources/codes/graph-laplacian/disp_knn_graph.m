function disp_knn_graph(I,J,z,f)

% display graph
hold on;
for t=1:length(I)
        j = J(t);
        i = I(t);
        c = (f(i)+f(j))/2;
        plot([z(i), z(j)], 'color', [c 0 1-c], 'LineWidth', 2);
end
for i=1:length(z)
    plot(z(i), '.', 'color', [f(i) 0 1-f(i)], 'MarkerSize', 20);
end

axis equal; 

end