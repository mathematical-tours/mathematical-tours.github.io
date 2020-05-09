function [fS,fH]=fusionca(C,S)
% fusionca	merge two increasing sequences.
%		fusionca(C,S) gives the resulting slopes fS and
%		indices fH where each slope support the epigraph.
%		All in all, it amounts to finding the first indice i
%		such that C(i-1)<=S(j)<C(i)
%
%		All slopes not needed to build the affine interpolation are
%		removed (SS may then be	smaller than S). For faster
%		computation use fusionma.m (it uses matlab built in
%		sort function to speed up computations).
%
% Usage:
%  [fS,fH]=fusionca(C,S)
%
% Input parameters:
%  C,S - Increasing column vectors. We want to know the conjugate at
%	S(j). C gives the slope of planar input points.
%
%
% Output parameter:
%  fS - Slopes supporting the epigraph. This is a subset of S. Column vector.
%  fH - Index at which the slope support the epigraph. It is a
%	two-column matrix. Indeed, a slope S(j) may touch at, at most, two
%	points (since C and S are increasing) which would be X(fH(j,1)) and
%	X(fH(j,1)). Of course, in most applications we only need to know one
%	point. However there are cases where one would like to know both points.
%
% Being called by:
%  fusion - Choose either of fusionca.m or fusionma.m
%
% Example:
%  X=[-5:0.5:5]'
%  Y=X.^2
%  C=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
%  S=C-0.25*ones(size(C))
%  [fS,fH]=fusionca(C,S)


% initialisation steps
n=size(C,1)-1;m=size(S,1);
i=1;j=1;
fS=[S(1);S(2)];

% We deal with the first slope
while S(j)>C(i), i=i+1; end

if S(j)==C(i)
  fH=[i i+1];
  i=i+1;
else
  fH=[i 0];
end

% We deal with the second slope
j=j+1;
while S(j)>C(i), i=i+1; end

if S(j)==C(i)
  fH(2,1:2)=[i i+1];
  i=i+1;
else
  fH(2,1:2)=[i 0];
end

j=j+1;
S(size(S,1)+1)=0; %Useful to make the last test of the loop

% Now we deal with the remaining slopes
while (j<=m & i<=n+1 & S(j)<=C(n+1)),
  while S(j)>C(i), i=i+1; end
  if fH(size(fH,1)-1,2)==0
	indice=fH(size(fH,1)-1,1);
  else
	indice=fH(size(fH,1)-1,2);
  end

  if S(j)==C(i)
	if indice<i
	  fH(size(fH,1)+1,:)=[i i+1];
	  fS(size(fS,1)+1)=S(j);
	  i=i+1;
	else
	  fH(size(fH,1),:)=[i i+1];
	  fS(size(fS,1))=S(j);
	  i=i+1;
	end
  elseif indice<i
	fH(size(fH,1)+1,:)=[i 0];
	fS(size(fS,1)+1)=S(j);
  else
	fH(size(fH,1),:)=[i 0];
	fS(size(fS,1))=S(j);
  end
j=j+1;
end

% Now we deal with the lastest slopes. Removing affine parts
% may allow us to skip some primal points.
if  (j<=m & i<=n+1 & S(j)>C(n+1))
  fH(size(fH,1)+1,:)=[size(C,1)+1 0];
  fS(size(fS,1)+1)=S(j);
  j=j+1;
end

if j<=m
  if fH(size(fH,1)-1,2)==0
	indice=fH(size(fH,1)-1,1);
  else
	indice=fH(size(fH,1)-1,2);
  end
  if indice<size(C,1)+1
	fH(size(fH,1)+1,:)=[size(C,1)+1 0];
	fS(size(fS,1)+1)=S(size(S,1)-1);
  else
	fH(size(fH,1),:)=[size(C,1)+1 0];
	fS(size(fS,1)+1)=S(size(S,1)-1);
  end
end

