N = 32;
x = [0:N/2, (N/2-1):-1:1]';
p = 9; i = 0;
for alpha = 0:pi/(2*(p-1)):pi/4
    i = i+1;
    y = fwt_interm(x,alpha);    
    subplot(1,5,i);
    plot(0:N-1,y,'k.:');
    axis square;
%    axis([0,N-1,-100,100]);
    axis tight;
    str = ['\alpha=', sprintf( '%.2f', alpha ) ];
    if(i==5) str = [str,' (Walsh)']; end;
    title( str );
end

saveas(gcf, '../transfo-walsh-interm-1', 'eps')
saveas(gcf, '../transfo-walsh-interm-1', 'png')
pause;

clf;
i = 0;
for alpha = (pi/4+pi/(2*(p-1))):pi/(2*(p-1)):pi/2
    i = i+1;
    y = fwt_interm(x,alpha);    
    subplot(1,5,i);
    plot(0:N-1,y,'k.:');
    axis square;
%    axis([0,N-1,-100,100]);
    axis tight;
    str = ['\alpha=', sprintf( '%.2f', alpha ) ];
    title( str );
end

saveas(gcf, '../transfo-walsh-interm-2', 'eps')
saveas(gcf, '../transfo-walsh-interm-2', 'png')