
addpath('toolbox/');
rep = 'results/';
[~,~] = mkdir(rep);

%%
% Read file.

filename =  'bovary.txt';
fid = fopen(filename, 'rt');
A = lower( fscanf(fid,'%c',Inf) );

%%
% Compute histogram.

S = 'abcdefghijklmnopqrstuvwxyz ''0123456789!()*,-.:;?';
x = zeros(size(A));
for i=1:length(S)
    h(i) = sum(A==S(i));
    x(A==S(i))=i;
end
x(x==0) = [];
h = h/sum(h);

clf; 
bar(1:length(S),h);
axis tight;
axis off;
set(gca, 'PlotBoxAspectRatio', [1 1/3 1]);
saveas(gcf, [rep 'hist.png'], 'png');

%%
% Compute Huffman tree.

T = compute_hufftree(h);
clf; plot_hufftree(T,0,S); 
set(gca, 'PlotBoxAspectRatio', [1 1/2 1]);
saveas(gcf, [rep 'tree.png'], 'png');

%%
% Encode the text.

Entropy = @(h)-sum(h.*log2(max(h,1e-20)));

[C,L] = huffman_gencode(T);
E = Entropy(h);
E1 = sum(h(:).*L);
% actual coding
b = perform_huffcoding(x(:),T,+1);
E2 = length(b)/length(x(:)); % should be equal to E1

fprintf(['Uniforme : ' num2str(ceil(log2(length(S)))) ' bits.\n']);
fprintf(['Entropie : ' num2str(E) ' bits.\n']);
fprintf(['Hufmann  : ' num2str(E1) ' bits.\n']);
