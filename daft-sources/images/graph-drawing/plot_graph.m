function plot_graph(A,xy)

% Dessine un graph en l'applattissant avec diverses
%   méthodes.
%
%   Copyright (c) 2003 Gabriel Peyré

global lang;
if ~strcmp(lang,'eng') % default is french
    lang = 'fr';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dessine de plusieurs façon un graph
n = length(A);
clf;
subplot(1,2,1);
gplot(A,xy,'k.-');
axis tight;
axis square;
axis off;
if strcmp(lang,'eng')==1
    title('"Manual" drawing');
else
    title('Dessin "Manuel"');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectral graph drawing : utilise les vecteurs propres
%   du laplacien.
lap = compute_laplacian(A);

disp('Performing SVD.');
tic;
[U,S,V] = svd(lap);
disp( sprintf('CPU time : %.2f.', toc) );
xy_spec = U(:,(n-2):(n-1));

% essaie de redresse le graphe
xy_spec = rectify_embedding(xy,xy_spec);

subplot(1,2,2);
gplot(A,xy_spec,'k.-');
axis tight;
axis square;
axis off;
if strcmp(lang,'eng')==1
    title('Combinatorial laplacian');
else
    title('Laplacien combinatoire');
end