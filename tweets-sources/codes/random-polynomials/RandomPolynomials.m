%%
% Plot zeros of random polynomials.

addpath('../toolbox/');
rep = MkResRep();

N = 200;
P = randn(N+1,1) + 1i * randn(N+1,1);
dg = (0:N)'; % degree

name = 'elliptic';
name = 'flat';
name = 'kac';

switch name
    case 'kac'
        % Kac polynomials in which {c_i := 1}.
        w = ones(N+1,1);
    case 'flat'
        % Flat polynomials (or Weyl polynomials) in which {c_i := \frac{1}{\sqrt{i!}}}.
        w = 1./sqrt(factorial(dg));
end

q = 70; 
nlist = round(linspace(3,N,q));
for it=1:q
    
    n = nlist(it);
    
    % Elliptic polynomials (or binomial polynomials) in which {c_i := \sqrt{\binom{n}{i}}}.
    if strcmp(name, 'elliptic')
        w = zeros(n+1,1);
        for i=0:n
            w(i+1) = sqrt(nchoosek(n,i));
        end
    end
    
    Q = w(1:n+1) .* P(1:n+1);
    z = roots(Q(end:-1:1));

    s = (it-1)/(q-1);
    clf; hold on;
    plot(z, '.', 'MarkerSize', 20, 'Color', [s 0 1-s]);
    plot(exp(2i*pi*linspace(0,1,100)), 'k', 'LineWidth', 1);
    axis equal; axis([-1 1 -1 1]*1.3); axis off;
    drawnow;
    
    % saveas(gca, [rep 'anim-' znum2str(it,2) '.png']);
end