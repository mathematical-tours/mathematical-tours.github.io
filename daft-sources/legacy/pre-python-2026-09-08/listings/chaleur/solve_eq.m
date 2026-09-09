function u = solve_eq(t, fcoef_val)
% précision du tracé
prec = 300; h = 1/prec;
u= zeros(prec+1, 1);
% taille de la solution
M = length(fcoef_val);
% fréquences des coefficients
v = [0:M/2,-M/2+1:-1]';
% calcule la solution
for i=0:prec
    x = i*h;
    w = exp(-2.0*pi*pi*t*v.*v + 2.0i*pi*x*v) .* fcoef_val;
    u(i+1) = sum(w);
end