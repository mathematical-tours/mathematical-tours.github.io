%%
% Test

X=[-5:0.5:5]'
Y=X.^2
S=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
[SS Conj]=LLTd(X,Y,S);

n = 256*4;
x = linspace(-1,1,n)';

q = 70;
for it=1:q
    u = (it-1)/(q-1);
    %
    f = x.^2 + .3*u*sin(6*pi*x);
    f = x.^2 + 2*u*x.^3;
    f = x.^2 - .3*u*abs(x-.5) - .2*u*abs(x+.5);
    s = linspace(-10,10,n*10)';
    [ss,fs]=LLTd(x,f,s);
    [x1,f1]=LLTd(ss,fs,x);
    %
    clf;
    subplot(2,1,1);
    plot(ss,fs);
    subplot(2,1,2);
    hold on;
    plot(x,f, 'b');
    plot(x1,f1, 'r');
    drawnow;
end