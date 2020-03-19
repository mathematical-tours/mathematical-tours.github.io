function [trset,teset ] = holdout( A,p )
n=size(A,1);
trdata=floor((n*p*0.01));
tedata=n-trdata;
index=randperm(n);
teset=zeros(tedata, size(A,2));
trset=zeros(trdata, size(A,2));
for i=1:tedata
    teset(i,:)=A(index(i),:);
end
for i=1:trdata
    trset(i,:)=A(index(i+tedata),:);
end



end

