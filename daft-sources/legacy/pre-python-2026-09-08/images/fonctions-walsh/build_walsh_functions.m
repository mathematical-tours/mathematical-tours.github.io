function build_walsh_functions(p)

global M;
global W;

N = 2^p;

W = ones(1,1);
for i=1:p
    W = [W,W;W,-W];
end

M = zeros(N,N);

for j=1:N
    t = 0;	% nombre de changements de signes
    s = W(1,j);
    for i=2:N
        if s~=W(i,j)
            t = t+1;
        end
        s= W(i,j);
    end
    M(:,t+1) = W(:,j);
end