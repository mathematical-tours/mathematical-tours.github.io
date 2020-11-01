[im, cm] = imread('maurice.png');
colormap(cm);

[N,P] = size(im);

U = double(im);

% constantes 
h = 0.0002;
d = 1/N;
s = h/d^2;
theta = 0.5;
nbr_row = 3;
nbr_col = 3;

display_row(U, s, 0, nbr_row, nbr_col, 1);
display_row(U, s, 0.5, nbr_row, nbr_col, 2);
display_row(U, s, 1, nbr_row, nbr_col, 3);


saveas(gcf, '../eq-chaleur-diff-finie', 'eps')
saveas(gcf, '../eq-chaleur-diff-finie', 'png')