%%
% Test of linear (and quadratic) discriminant analysis.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);

k = 8; % number of classes
n = 200; % point per classes

% initial positions/rotation/scale/aniso
m = rand(k,1)+1i*rand(k,1);
r = rand(k,1);
s = ( .8+.2*rand(k,1) ) * .1;
a = .25 * ones(k,1); % .25+.25*rand(k,1);

lda = 1;
r = r*0+.3;
if not(lda)
    r = rand(k,1);
end
a = a*0 + 1/4;
s = s*0 + .12;

r = r*0-.3;
a = a*0 + 1/8;
s = s*0 + .1;

x0 = randn(k,n)+1i*randn(k,n);
x = exp(2i*pi*r) .* ( real(x0).*s + 1i*imag(x0).*a.*s  ) + m;

co = distinguishable_colors(k+1); co(3,:) = [];
clf; hold on;
for i=1:k
    plot( real(x(i,:)), imag(x(i,:)), '.', 'color', co(i,:));
end
l = .2;
axis equal;
axis([-l,1+l,-l,1+l]);


% animate points in a square with bouncing
rand('state', 123); randn('state', 123);

m = rand(k,1)+1i*rand(k,1);
eta = .02; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);

% grid points
l = .2;
n0 = 256;
t0 = linspace(-l,1+l,n0);
[Y,X] = meshgrid(t0,t0);
Z = X+1i*Y;
Z = transpose(Z(:));
co = [[1 0 0];[0 1 0];[0 0 1];[0 0 0];[1 0 1];[0 1 1];[1 1 0];[1/2 1/2 1]];

q = 120;
for it=1:q
    t = (it-1)/(q-1);
    rt = r;
    if ~lda
        rt = r + (-1).^(1:k)' * t/2;
    end
    % position points
    x = exp(2i*pi*rt) .* ( real(x0).*s + 1i*imag(x0).*a.*s  ) + m;    
    % classification, via max of probability
    U = exp(-2i*pi*rt) .* (Z-m);
    U = real(U)./s + 1i*imag(U)./(s.*a);
    % S = 1./(s.^2.*a) .* exp( -abs(U).^2/2 );
    S = -log(s.^2.*a) - abs(U).^2/2;
    [~,I] = max(S, [], 1);
	% turn into color array 
    C = zeros(n0*n0,3);
    for i=1:k
        for j=1:3
            C(I==i,j) = .6+.6*co(i,j);
        end
    end
    C = reshape(C,n0,n0,3);
    %
    clf; hold on;
    imagesc(t0,t0,permute(C,[2 1 3]));
    for i=1:k
        plot( real(x(i,:)), imag(x(i,:)), '.', 'color', .7*co(i,:));
    end
    axis equal; axis off;
    axis([-l,1+l,-l,1+l]);
    drawnow;
    mysaveas(it);
    
    % DO HERE STUFF
    m = m + v;
    I = find( real(m)<0 | real(m)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(m)<0 | imag(m)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));
end
