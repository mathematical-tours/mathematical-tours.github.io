function [bbx,bby] = bb (X,Y)
% bb		Compute the planar convex hull of the set (X(i),Y(i)).
%		By requesting X to be sorted we obtain a linear-time
%		algorithm.
%
% Usage:
%  [bbx,bby] = bb (X,Y)
%
% Input parameters:
%  X,Y - Column vectors of same size. (X(i),Y(i)) is the planar set
%	whose convex hull we compute. X must be sorted increasingly
%	with distinct points: X(i)<X(i+1).
%
%
% Output parameter:
%  bbx,bby - Verteces of the convex hull. Obviously bbx is a subset of
%	X and bby a subset of Y. Column vectors.
%
%
% Being called by:
%  LLTd - Compute the discrete Legendre transform.
%
% Example:
%  X=[-2:0.25:2]'
%  Y=(X.^2-ones(size(X))).^2
%  [bbx,bby]=bb(X,Y)
%  plot(X,Y,'x');hold on;plot(bbx,bby,'o');
%
% bb uses the Beneat-Beyond algorithm [Preparata and Shamos,
% Computational Geometry, Springer Verlag, 1990]

% Initialisation
n=size(X,1);
I=[1:n];v=2;CX=NaN*ones(n,1);CY=CX;
CX(1:2)=[X(1);X(2)];CY(1:2)=[Y(1);Y(2)];
v=2;

% We look at each point only once, this gives the linear-time
% algorithm.
for i=3:n,
  y=((CY(v)-CY(v-1))/(CX(v)-CX(v-1)))*(X(i)-CX(v))+CY(v);
  while (v>2 & Y(i)<=y),
	% We erase points which are not vertices of the convex hull
	v=v-1;
	y=((CY(v)-CY(v-1))/(CX(v)-CX(v-1)))*(X(i)-CX(v))+CY(v);
  end
  if v>2
	if Y(i)==y
	  CX(v)=X(i);CY(v)=Y(i);
	else
	  CX(v+1)=X(i);
	  CY(v+1)=Y(i);
	  v=v+1;
	end
  else
	if Y(i)>y % Trivial convex hull
	  CX(3)=X(i);CY(3)=Y(i);v=3;
	else
	  CX(2)=X(i);CY(2)=Y(i);
	end
  end;
end
bbx=CX(1:v);bby=CY(1:v);

