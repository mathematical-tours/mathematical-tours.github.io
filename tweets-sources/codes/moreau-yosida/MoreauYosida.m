rep = '../results/moreau-yosida/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);



name = 'cvx';
name = 'ncvx';
name = 'abs';

switch name
    case 'abs'
        f = @(x)abs(x);
    case 'cvx'
        f = @(x)abs(x-.3) + 2*abs(x+.3);
    case 'ncvx'
        f = @(x)(x<=0).*( .3-abs(x+.3) ) + (x>=0).*( -.3+abs(x-.3) );
end


n = 1001;
p = 1024*50;
x = linspace(-1,1,n);
y = 3*linspace(-1,1,n);

% f_eps(x) = umin_{y} f(y) + 1/(2*epsilon)*|x-y|^2
[Y,X] = meshgrid(y,x);
MR = @(f,epsilon,x)min(f(Y)+1/(2*epsilon)*abs(X-Y).^2,[],2);

q = 10;
tlist = linspace(.5,0,q);
clf; hold on;
for i=1:length(tlist)
    t = tlist(i);
    c = (i-1)/(length(tlist)-1);
    if t==0
        u = f(x);
    else
        u = MR(f,t);
    end
    plot(x, u, 'LineWidth', 2, 'Color', [c 0 1-c]);
end
axis tight; box on;
SetAR(1/2); 
set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, [rep name '-regul.eps'], 'epsc');
