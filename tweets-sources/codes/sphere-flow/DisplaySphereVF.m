function DisplaySphereVF(Vh,Fh, V,g)

% display the vector field 
rho = 1.01; Vr = rho*V; % base point for drawing;
xi = .07;
hold on;
plot_mesh(Vh,Fh);
for i=1:size(Vr,2)
    m = [Vr(:,i), Vr(:,i) + xi*g(:,i)];
    plot3(m(1,:), m(2,:),m(3,:),'k', 'LineWidth', 2); 
end
% plot3(Vr(1,:), Vr(2,:), Vr(3,:), 'k.', 'MarkerSize', 8);


end
