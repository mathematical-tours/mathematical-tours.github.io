% test du décodage de R(1,n)
n = 6;
d = 2^n;                        % taille
e = 2^(n-2);                    % nbr d'erreurs max
a = floor( rand*d );            % numéro à encoder
b = rand>0.5;                   % signe

f = encode_rm(a,b,n);
err = floor( 1+rand(e,1)*d );   % e positions tirées au hasard
g = f; g(err) = 1-g(err);       % vecteur perturbé
gg = decode_rm(g);
norm( f-gg )                    % doit valoir zéro