function infconv = infconvd(Xf,Yf,Xg,Yg,S)
% infconvd	compute the inf-convolution of two convex functions f
%		and g satisfying Yf=f(Xf) and Yg=g(Xg). The
%		inf-convolution is computed by using the formula:
%		infd(f,g)(x):=inf_{y}[f(y)+g(x-y)]=(f^*+g^*)^*.
%
%
% Usage:
%  infconv = infconvd(Xf,Yf,Xg,Yg,S)
%
% Input parameters:
%  Xf,Yf - Column vectors. Usually Yf(i)=f(Xf(i))
%  Xg,Yg - Column vectors. Usually Yg(i)=g(Xg(i))
%  S - Column vector. We want to know the infconvolution at S(j)
%
%
% Output parameter:
%  infconv - Infconvolution of both functions. This is a column
%	vector.
%
% Functions called:
%   LLTd - Compute the discrete Legendre transform.
%
%
% Example:
%  X=[-5:0.5:5]';S=X;
%  Yf=0.5*X.^2;Yg=Yf;
%  infconv=infconvd(X,Yf,X,Yg,S)

m=size(S,1);n=size(Xf,1);
[Sf,Cf]=LLTd(Xf,Yf,S);
[Sg,Cg]=LLTd(Xg,Yg,S);

% Check that Sf=Sg
if (size(Sf,1)~=size(S,1))|(size(Sg,1)~=size(S,1)),
  error('Error : use fusionma.m instead of fusionca.m in fusion.m');
end;

X=sort([Xf;Xg]);
Y=Cf+Cg;
% The sort function is faster than a
% for-loop. Anyway, by using a for-loop,
% sorting can be done in linear-time.
[SS,infconv]=LLTd(S,Y,X);

