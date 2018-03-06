function g = apply_colormap(f,CM)
% g = apply_colormap(f,CM);

N = size(CM,1);
f = floor(rescale(f)*(N-1))+1;
g = reshape(CM(f,:), [size(f) 3]);

end