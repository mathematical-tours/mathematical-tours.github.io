function w = walsh(n)

global M;
global N;

% calcul de (i,j)
i = 1;
j = 1;
s = 1;
while s<n
    s = s+1;
    j = j+1;
    i = i-1;
    if j==N+1
        j = i+2;
        i = N;
    end
    if i==0
        i=j;
        j=1;
    end
end

w = M(:,j)*M(:,i)';