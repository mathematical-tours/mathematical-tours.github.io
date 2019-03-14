function mybox(U, col, lw)


x = U(1:2); 
y = U(3:4);

plot( [x(1)+1i*y(1),x(2)+1i*y(1),x(2)+1i*y(2),x(1)+1i*y(2),x(1)+1i*y(1)],'LineWidth',lw, 'Color', col);

end