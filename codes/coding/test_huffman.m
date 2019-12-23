%%
% test for huffman coding of images.

addpath('toolbox/');

rep = 'results/huff/';
if not(exist(rep))
    mkdir(rep);
end

%%
% Global parameters.

% name of the image
name = 'boat'; 
name = 'hibiscus';
% #quantification
if not(exist('q'))
    q = 4;
end
qmax = q;
fid = fopen([rep 'q' num2str(q) '-results.txt'], 'wt');
% image size 
n = 256;

%%
% Helpers.

normalize = @(x)x/sum(x(:));
myhist = @(x,y)normalize(hist(x,y));
Quant = @(x,q)min(floor( rescale(x)*q  ), q-1);
Entropy = @(h)-sum(h.*log2(max(h,1e-20)));
myaxis = @(x,h)axis([min(x)-1/2 max(x)+1/2 0 max(h)*1.02]);
SetAR = @()set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'FontSize', 20);

%% 
% Load image.

x0 = load_image(name, n);
x0 = rescale(sum(x0,3));
imwrite(x0, [rep 'original.png']);

%%
% Quantization.

x = Quant(x0,q);
imwrite(rescale(x), [rep 'q' num2str(q) '-image.png']);

%%
% Coding pixels in {0,...,q-1}. 

hx = myhist(x(:),0:q-1);
Tx = compute_hufftree(hx);

clf;
imageplot({x0 x});

clf; 
bar(0:q-1, hx); axis tight;
myaxis(0:qmax-1, hx);
SetAR();
saveas(gcf, [rep 'q' num2str(q) '-pxl-histo.eps'], 'epsc');

clf; plot_hufftree(Tx, -1); 
saveas(gcf, [rep 'q' num2str(q) '-pxl-tree.eps'], 'epsc');

[Cx,Lx] = huffman_gencode(Tx);
bx = perform_huffcoding(x(:)+1,Tx,+1);
Ex = Entropy(hx);
Ex1 = sum(hx(:).*Lx);
Ex2 = length(bx)/length(x(:)); % should be equal to Ex1
fprintf('Entropy=%.3f, Huffman=%.3f, log2(q)=%.3f\n', Ex, Ex1, log2(q));

fprintf(fid, '----- PIXELS ------\n');
fprintf(fid, 'Entropy=%.3f, Huffman=%.3f, log2(q)=%.3f\n', Ex, Ex1, log2(q));
for i=1:length(Cx)
    fprintf(fid, '%d ---> %s\n', i, num2str(Cx{i}));
end

% Store in a file.
fidC = fopen([rep 'q' num2str(q) '-plx-code.txt'], 'wt');
fprintf(fidC, '%d',bx(:)); 
fclose(fidC);

%%
% Differential coding in {0,...,2*q-2}.

y = diff(x(:))+q-1;
hy = myhist(y(:),0:2*q-2);

clf; 
bar(-(q-1):q-1, hy); axis tight;
myaxis(-(qmax-1):qmax-1, hy);
SetAR();
saveas(gcf, [rep 'q' num2str(q) '-diff-histo.eps'], 'epsc');

Ty = compute_hufftree(hy);
[Cy,Ly] = huffman_gencode(Ty);
by = perform_huffcoding(y(:)+1,Ty,+1);

clf; plot_hufftree(Ty, -q); 
saveas(gcf, [rep 'q' num2str(q) '-diff-tree.eps'], 'epsc');

Ey = Entropy(hy);
Ey1 = sum(hy(:).*Ly);
Ey2 = length(by)/length(y); % should be equal to Ey1

fprintf('Entropy=%.3f, Huffman=%.3f, log2(q)=%.3f\n', Ey, Ey1, log2(q));

fprintf(fid, '----- DIFFERENCES ------\n');
fprintf(fid, 'Entropy=%.3f, Huffman=%.3f, log2(q)=%.3f\n', Ey, Ey1, log2(q));
for i=1:length(Cy)
    fprintf(fid, '%d ---> %s\n', i-q, num2str(Cy{i}));
end
fclose(fid);

% Store in a file.
fidC = fopen([rep 'q' num2str(q) '-diff-code.txt'], 'wt');
fprintf(fidC, '%d',bx(:)); 
fclose(fidC);