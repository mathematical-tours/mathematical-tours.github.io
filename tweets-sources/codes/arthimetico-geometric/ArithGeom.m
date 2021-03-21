addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);


n = 256*2;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);

clf; mydisp((X+Y)/2);
mysaveas('arith', 1);
clf; mydisp(sqrt(X.*Y));
mysaveas('geom', 1);


A = X; G = Y;
for it=1:100
    [A,G] = deal( (A+G)/2,  sqrt(A.*G) );
end
M = A;

clf; mydisp(M);
mysaveas('arith-geom', 1);