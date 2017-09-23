%%
% Simple test.

addpath('toolbox/');

h = [5 5 2 2 3 1];
h = h/sum(h);
T = compute_hufftree(h);
clf; plot_hufftree(T); 

[C,L] = huffman_gencode(T);
for i=1:length(C)
    fprintf('%d --> %s\n', i, num2str(C{i}));
end

b  = perform_huffcoding([1 2 1 3],T,+1);
perform_huffcoding(b,T,-1)