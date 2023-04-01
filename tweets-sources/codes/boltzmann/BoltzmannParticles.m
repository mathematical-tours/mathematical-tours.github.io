%%
% Evolves particles according to interaction laws.

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


xmax = 1/2;

kappa = .04; % collision radius


k = 30;
k = 500;
k = 100;
kappa = .4 * xmax / sqrt(k);


S = 1000;
X = poissonDisc([1 1]*xmax*S,1.5 * S * kappa,k)/S;
x = X(:,1) + 1i*X(:,2);
if length(x)<k
    error('not able to seed particles');
end

eta = .001/5; v = randn(k,1) + 1i*randn(k,1); v = eta*v./abs(v);

% display
sz = kappa * 2500; % 75; % size
sz = kappa * 2000; % 75; % size
c = rand(k,3); % color




dotp = @(x,y)real(x*conj(y));

q = 120*200;
jt = 0;
for it=1:q
    % DO HERE STUFF
    x = x + v;
    I = find( real(x)<0 | real(x)>1 );
    v(I) = -real(v(I)) + 1i*imag(v(I));
    I = find(  imag(x)<0 | imag(x)>1 );
    v(I) = real(v(I)) - 1i*imag(v(I));


    D = abs(x-transpose(x));
    [B,A] = meshgrid(1:k,1:k);
    I = find(D<kappa & B>A);
    A = A(I); B = B(I);
    % treat each collision 
    for num=1:length(B)
        a = A(num); b = B(num);
        u = (v(a)-v(b))/abs(v(a)-v(b)); u1 = 1i*u; % tangent/orthogonal to collision direction
        v([a,b]) = [ dotp(v(a),u1)*u1 + dotp(v(b),u)*u; ...
                     dotp(v(b),u1)*u1 + dotp(v(a),u)*u ]; 
    end

    if mod(it,40)==1
        clf;    
        clf; scatter( real(x), imag(x),sz, c, 'filled' );
        axis equal;
        e = kappa/2;
        axis([-e 1+e -e 1+e]); axis on; box on;
        set(gca, 'Xtick', [], 'Ytick', []);
        drawnow;
        jt = jt+1;
        mysaveas(jt);
    end
end




% AutoCrop(rep, 'simul-');