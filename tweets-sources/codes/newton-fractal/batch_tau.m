rho_list = 0:.05:1.5;

for irho=1:length(rho_list)
    rho = rho_list(irho);
    newton_fractal;
end