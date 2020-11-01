clear;
clf;
global N;
global M;
global W;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dessine les matrices des fonctions 1D ordonnées de 2 façons

p = 6;
N = 2^p;
build_walsh_functions(p);

subplot(1,2,1);
image((W+1)*128);
title('Ordre naturel');
colormap(gray);
axis image;
grid on;
set(gca,'GridLineStyle', '-');
set(gca, 'XTick', (1:N)+0.5);
set(gca, 'YTick', (1:N)+0.5);
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');

subplot(1,2,2);
image((M+1)*128);
title('Ordre par nombre de changements de signes');
colormap(gray);
axis image;
grid on;
set(gca,'GridLineStyle', '-');
set(gca, 'XTick', (1:N)+0.5);
set(gca, 'YTick', (1:N)+0.5);
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');

saveas(gcf, '../fonctions-walsh', 'eps');
saveas(gcf, '../fonctions-walsh', 'png');

pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dessine une compression progressive d'une fonction 1D
clf;

p = 8;
N = 2^p;
build_walsh_functions(p);


x = (0:N-1)/N;
f = cos(1*exp(4.5*x')) .* cos((x'-1/2)*pi).^2;

compression = [100,80,60,40,30];
for i=1:5
    subplot(1,5,i);
    plot( 0:N-1, decompose(f, compression(i)/100*N) );
    axis([0,N-1, -1,1]);
    axis square;
    str = sprintf('%d%%', compression(i) );
    title(str);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);        
end
saveas(gcf, '../compression-walsh-1d-1', 'eps');
saveas(gcf, '../compression-walsh-1d-1', 'png');

compression = [25,20,15,10,5];
for i=1:5
    subplot(1,5,i);
    plot( 0:N-1, decompose(f, compression(i)/100*N) );
    axis([0,N-1, -1,1]);
    axis square;
    str = sprintf('%d%%', compression(i) );
    title(str);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);        
end
saveas(gcf, '../compression-walsh-1d-2', 'eps');
saveas(gcf, '../compression-walsh-1d-2', 'png');

pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dessine les différents spectres obtenus.
s1 = W*f;
s2 = M*f;
subplot(2,1,1);
plot(0:N-1, abs(s1), 'k');
title('Spectre, ordre naturel');
axis([0,N-1,0,max(s1)]);
subplot(2,1,2);
plot(0:N-1, abs(s2), 'k');
title('Spectre, ordre par nombre de changements de signe');
axis([0,N-1,0,max(s2)]);

saveas(gcf, '../spectre-walsh-1d', 'eps');
saveas(gcf, '../spectre-walsh-1d', 'png');

pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classement des fonctions de Walsh 2D.

clf;
p = 3;
N = 2^p;
build_walsh_functions(p);

for i=1:N^2
    subplot(N,N,i);
    image((walsh2d(i)+1)*128);
    colormap(gray);
    axis image;
    grid on;
    set(gca,'GridLineStyle', '-');
    set(gca, 'XTick', (1:N)+0.5);
    set(gca, 'YTick', (1:N)+0.5);
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
end

saveas(gcf, '../fonctions-walsh-2d', 'eps');
saveas(gcf, '../fonctions-walsh-2d', 'png');

pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dessine une compression progressive d'une image 2D
p = 6;
N = 2^p;
W = ones(1,1);
build_walsh_functions(p);

[f, cm] = imread('image.bmp');
f = double(f);

compression = [100,80,60,50,40];
for i=1:length(compression)
    subplot(1,5,i);
    image( decompose2d(f, compression(i)/100*N^2)*255 );
    str = sprintf('%d%%', compression(i) );
    title(str);
    axis image;
    grid on;
    colormap gray;
   
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end
saveas(gcf, '../compression-walsh-2d-1', 'eps');
saveas(gcf, '../compression-walsh-2d-1', 'png');

compression = [30,20,15,10,5];
for i=1:length(compression)
    subplot(1,5,i);
    image( decompose2d(f, compression(i)/100*N^2)*255 );
    str = sprintf('%d%%', compression(i) );
    title(str);
    axis image;
    grid on;
    colormap gray;
   
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
end
saveas(gcf, '../compression-walsh-2d-2', 'eps');
saveas(gcf, '../compression-walsh-2d-2', 'png');

pause;