n = 500*2;
a0 = rand(n,1); 
a1 = rand(n,1); 

x = linspace(0,1,n)'; 
a0 = (abs(x-.3)<.1) + 1e-10;
a1 = (abs(x-.7)<.17) + 1e-10;

a0 = a0/sum(a0);
a1 = a1/sum(a1);


q = 50;
for it=1:q
    t = (it-1)/(q-1);
    rev = @(x)interp1(x,(0:n)/n,(0:n)/n,'linear','extrap'); % inversion
    at = diff( rev( (1-t)*rev([0;cumsum(a0)]) + t*rev([0;cumsum(a1)]) ) );
    clf; hold on;
    plot(x, [a0 a1]);
    plot(x, at, 'k');
    drawnow;
end