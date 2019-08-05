%%
% Plot Gaussian scale mixtures.


addpath('../toolbox/');
rep = MkResRep();


A = 7; 
n = 1001;
t = linspace(-A,A,n)';

G = @(z)1/sqrt(2*pi*z) * exp(-t.^2/(2*z));

r = 200;
Z = linspace(0.01,8,r)';

nu = 2; 

% scales
sc = Z.^(-nu/2-1) .* exp(-nu./(2*Z));
sc = sc / sum(sc);
clf; plot(Z, sc, '.-');

Z = 1./gamrnd(gamma*ones(r,1),gamma*ones(r,1));
sc = ones(r,1)/r;

clf; hold on;
H = zeros(n,1);
for iz=1:r
    s = (iz-1)/(r-1);
    z = Z(iz);
    h = G(z) * sc(iz);
    H = H + h;
    plot( h, 'Color', [s 0 1-s] );
end

normalize = @(H)H / (sum(H)*2*A/n);

H = H / (sum(H)*2*A/n);
% student
T = (1+t.^2/nu).^(-(nu+1)/2);
T = T / (sum(T)*2*A/n);

clf; hold on;
plot( H, 'b' );
plot( T, 'k--' );

return;

r = 6;
nu_list = linspace(.1,6,r);
nu_list = .05 + 4*linspace(0,1,r).^2;
B = 4;
Z = linspace(1e-3,B,n);

clf; hold on;
for it=1:r
    s = (it-1)/(r-1);
    nu = nu_list(it);
    h = Z.^(-nu/2-1) .* exp(-nu./(2*Z));
    h = h / (sum(h)*B/n);
    plot( Z, h, 'Color', [s 0 1-s], 'LineWidth', 2 );
    box on;
end
set(gca, 'Fontsize', 20);
saveas(gcf, [rep 'scale.eps'], 'epsc');

clf; hold on;
for it=1:r
    s = (it-1)/(r-1);
    nu = nu_list(it);
    h = (1+t.^2/nu).^(-(nu+1)/2);
    h = h / (sum(h)*2*A/n);
    plot( t, h, 'Color', [s 0 1-s], 'LineWidth', 2 );
    axis tight;
    box on;
end
set(gca, 'Fontsize', 20);
saveas(gcf, [rep 'pdf.eps'], 'epsc');


