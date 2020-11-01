clear;
clf;

subplot(1,4,1);
draw_matrix(5);
subplot(1,4,2);
draw_matrix(6);
subplot(1,4,3);
draw_matrix(7);
subplot(1,4,4);
draw_matrix(8);


saveas(gcf, '../rev-bit-matrix', 'eps');
saveas(gcf, '../rev-bit-matrix', 'png');