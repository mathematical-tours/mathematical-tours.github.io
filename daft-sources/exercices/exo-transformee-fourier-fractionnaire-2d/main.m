% exemples de transformée de Fourier fractionnaire (1D).

n = 256;
h = -1:2/(n-1):1;


[Y,X] = meshgrid( (0:(n-1))*2, (0:(n-1))*2 );
x = sin((X+Y));

[Y,X] = meshgrid( h, h );
x = (X.^2+Y.^2)<0.1^2;

subplot( 3,3,1 );  
x1 = (X.^2+Y.^2)<0.3^2;
imagesc(-x1);
title('Image d''origine')
axis image
axis off;
colormap gray;




k=1;

for beta = 0.6:0.2:2.0 % 0.4:0.3:2.8

    k = k+1;
    
    y = frft2d( x, beta );
    % y = fftshift(y);
    y = abs(y);
    m1 = min(min(y));
    m2 = max(max(y));
    sc = 10000; %100000;
    y = sc*(y-m1)/(m2-m1);
    y = 256-y;
    
    subplot( 3,3,k );
    
    image(y);
    axis image
    axis off;
    colormap gray;
    
    title(sprintf('\\alpha=%.2f', beta));

end

saveas(gcf, '../transformee-fourier-fractionnaire-2d', 'eps')
saveas(gcf, '../transformee-fourier-fractionnaire-2d', 'png')