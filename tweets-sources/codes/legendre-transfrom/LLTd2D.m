function Conj2D=LLTd2D(X1,X2,Y,S1,S2)
% LLTd2D	Compute numerically the discrete Legendre transform of
%		a function defined on the plane: Y=f(X1(i),X2(j))
%		at slopes (S1(l),S2(k)), i.e., it gives
%		Conj2D(l,k)=max_{i,j}[S1(l)*X1(i)+S2(k)*X2(j)-Y(i,j)].
%		It uses the univariate Linear-time Legendre Transform
%		algorithm with the factorization
%		f^*(s1,s2)=sup_x1[s1 x1 + sup_x2 [s2 x2 -f(x1,x2)]]
%		with Y(i1,i2)=f(X1(i1),X2(i2)), to speed up
%		computations.
%
% Usage:
%  Conj2D=LLTd2D(X1,X2,Y,S1,S2)
%
% Input parameters:
%  X1,X2 - Column vectors. Usually Y(i1,i2)=f(X1(i1),X2(i2)).
%  Y - Matrix of real numbers.
%  S1,S2 - Column vectors. We want to know the conjugate at (S1(l),S2(k)).
%
%
% Output parameter:
%  Conj2D - numerical values of the conjugate at (S1(l),S2(k)). This is
%	a matrix.
%
% Functions called:
%  LLTd - Compute the univariate discrete Legendre transform.
%
% Example:
%  X1=[-5:0.5:5]'
%  X2=X1
%  Y=(X1*X2').^2
%  S1=[-200:20:200]';
%  S2=S1
%  Conj2D=LLTd2D(X1,X2,Y,S1,S2)
%  mesh(S1,S2,Conj2D)
%

n1=size(X1,1);n2=size(X2,1);
m1=size(S1,1);m2=size(S2,1);

% loop for each x1 in X1
for i1=1:n1,
  [S,Vt]=LLTd(X2,Y(i1,:),S2);
  Vs2(i1,:)=-Vt';
end;

% loop for each s2 in S2
for j2=1:m2,
  [S,Vt]=LLTd(X1,Vs2(:,j2),S1);
  Conj2D(:,j2)=Vt;
end;
