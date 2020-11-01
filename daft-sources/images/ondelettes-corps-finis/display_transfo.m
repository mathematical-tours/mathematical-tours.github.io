function display_transfo(f,psi, nbrow, row)

p = length(psi);
t = transfo_ondelettes(f,psi);

subplot(nbrow, 3, 1+3*(row-1));
plot(f, 'k.:');
axis square; axis tight;
if(row==1) title('Fonction f'); end;

subplot(nbrow, 3, 2+3*(row-1));
p0 = (p-1)/2;
x = -p0:p0; 
c = [(p0+2):p, 1:(p0+1)];
plot(x, psi(c), 'k.:');
axis square; axis tight;
if(row==1) title('Fonction \psi'); end;

subplot(nbrow, 3, 3+3*(row-1));
imagesc(-abs(t));
colormap gray;
axis image;
xlabel('b');
ylabel('a');
if(row==1) title('Transformée W(f)(b,a)'); end;
% axis off;