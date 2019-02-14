L = .27;
k = 1-3/20;

k = 1-1/20;


addpath('../toolbox/');
rep = MkResRep();


t = 2*pi*linspace(0,25,2048*2);



Spiro = @(k,L)(1-k)*cos(t) + L*k*cos((1-k)/k*t) + 1i * ( (1-k)*sin(t) - L*k*sin((1-k)/k*t) );

Q = 70;

% k=r/R<1
klist = 1-19/25 + 0*linspace(0.5,1,Q); % 6/25
% d/R=L*k --> d=L*r>1
Llist = linspace(0,10,Q);


NN = @(x)x/max(abs(x(:)));

for i=1:Q
    c = (i-1)/(Q-1);
    clf;
    plot(NN(Spiro(klist(i),Llist(i))), 'LineWidth', 2, 'Color', [c 0 1-c]);
    axis([-1 1 -1 1]);
    axis off;
    drawnow;
    saveas(gcf, [rep 'spiro-' znum2str(i,2) '.png'], 'png');
end