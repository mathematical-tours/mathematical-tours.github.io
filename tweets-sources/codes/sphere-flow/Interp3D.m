function g = Interp3D(f,V)

n = size(f,1);

t = linspace(-1.1,1.1,n); g = [];
for i=1:3
    g(i,:) = interp3(t,t,t,f(:,:,:,i), V(1,:), V(2,:), V(3,:));
end
g = g - sum(g.*V) .* V;  % project on the tangent space sphere
g = g ./ sqrt(sum(g.^2)); 

end
