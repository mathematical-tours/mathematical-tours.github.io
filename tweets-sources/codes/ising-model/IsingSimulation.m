%%
% Test of ising model.
% for J>0, Ising model of a ferromagnet
% https://fr.mathworks.com/matlabcentral/fileexchange/62194-ising-model-and-metropolis-algorithm
% E = -J sum_ij sigma_i sigma_j


addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

n = 32*2;

J = 1;

Tc = 2*J / log(1+sqrt(2)); % Curie temperature, phase transition of the ferromagnet

T = Tc*2; % temperature

it = 0;
for a=1:10
    u = false;
    while true
        S = metropolis( sign(randn(n)) , T, J);
        if mean(S(:))>-.9 && mean(S(:))<.9
            break; % avoid pure black or white
        end
    end
    for b=1:5
        it = it+1;
        clf;
        imagesc(S); axis off; axis image;
        colormap(gray(256));
        drawnow;
        mysaveas(it);
    end
end
