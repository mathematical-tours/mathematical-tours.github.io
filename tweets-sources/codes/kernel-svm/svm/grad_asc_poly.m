function [alpha] = grad_asc_poly( X_Data,Y_Data,C,deg )

N=size(X_Data,1);
alpha=zeros(N,1);
alp2=zeros(N,1);
norm1=10e2;
tol=10e-3;
Ker=(X_Data*X_Data').^deg;

while norm1>tol
alp_old=alpha;
w1=(alp_old.*Y_Data).*X_Data;

for i=1:N
    eta=1/Ker(i,i);
    for j=1:N
        alp2(j,1)=alpha(j,1)*Y_Data(j,1)*Ker(i,j);
    end
    alpha(i,1)= alpha(i,1)+eta*(1-(Y_Data(i,1)*sum(alp2)));
    
    if alpha(i,1)<0
        alpha(i,1)=0;
    elseif alpha(i,1)>C
        alpha(i,1)=C; 
    end
    
end

    w2=(alpha.*Y_Data).*X_Data;
    norm1=norm(w2-w1);
end

end

