n = 16/2;
nmax = 1;
H = zeros(n);
for i=1:n
    x = zeros(n,1); x(i)=1;
    H(:,i) = Haar(x,nmax);
end

imagesc(H);
colormap gray;
axis image; axis off;

U = zeros(2*n);
for i=1:n
    U(i,2*i-1:2*i) = [1 1];
    U(n+i,2*i-1:2*i) = [1 -1];    
end
imagesc(U);
colormap gray;
axis image; axis off;