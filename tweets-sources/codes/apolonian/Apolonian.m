
addpath('../toolbox/');

name = 'elephant';
name = 'cat';
name = 'bunny';
name = 'france';

shape = '1';
shape = 'ncv';
shape = 'inf';
shape = '2';

randrot = 0;

rep = MkResRep(name);

switch shape
    case '2'
        Shp = @(u)abs(u);
    case '1'
        Shp = @(u)abs(real(u))+abs(imag(u));
    case 'inf'
        Shp = @(u)max(abs(real(u)),abs(imag(u)));
    case 'ncv'
        a = .4;
        Shp = @(u)( abs(real(u)).^a+abs(imag(u)).^a ).^(1/a);
end


n = 1000;
W0 = double(rescale(sum(load_image(name,n),3))>.5);
[Y,X] = meshgrid(1:n,1:n);
P = X(:)+1i*Y(:);

i = [2:n n]; j = [n 1:n-1];
cvl = @(F)( F+F(i,:)+F(j,:)+F(:,i)+F(:,j)+F(i,i)+F(j,j)+F(i,j)+F(j,i) )/9;
bnd = @(F)find( cvl(F)>0 & cvl(F)<1 & F==1 );

Dmat = @(x,y)Shp( x(:)*ones(1,length(y)) - ones(length(x),1)*transpose(y(:)) );

q = 80; % #shape

s = 10*1e3; % #random sample

rrand = @(m,v)floor(rand(m,1)*v)+1;

Col = distinguishable_colors(q+10);
[~,I] = sort(sum(Col,2), 'descend');
Col = Col(I(1:q),:);

m = 100; % #random samples
W = W0;
% for rendering
F = repmat(W0,[1 1 3]);
F = ones(n,n,3);
for it=1:q
    % random rotation factor
    om = 1;
    if randrot
        om = exp(2i*pi*rand);
    end
    % compute ext boundary
    B = bnd(W);
    % distance from inside to bound
    I = find(W==0); % inside
    I = I( rrand(min(length(I),s),length(I)) );
    %
    D = Dmat( om*P(B),om*P(I) );
    [d,K] = min(D,[],1);
    [d0,i] = max(d); i = I(i);
    % add the disk
    J = find(Shp(om*(P(i)-P))<d0);
    W(J) = 1;
    % rendering
    F(J) = Col(it,1);
    F(J+n^2) = Col(it,2);
    F(J+2*n^2) = Col(it,3);
    %
    imageplot(F);
    drawnow;
    % 
    imwrite(F, [rep name '-' shape '-' znum2str(it,2) '.png']);
end