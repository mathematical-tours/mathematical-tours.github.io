function [fS,fH]=fusionma(C,S)
% fusionma	merge two increasing sequences.
%		fusionma(C,S) gives the resulting slopes fS and
%		indices fH where each slope support the epigraph.
%		All in all, it amounts to finding the first indice i
%		such that C(i-1)<=S(j)<C(i)
%
% 		This function uses matlab syntax to speed up
%		computation, see fusionca.m for a more classical
%		programming of the same function
%
%		All slopes not needed to build the affine interpolation may be
%		removed by using fusionca.m instead of fusionma.m but
%		computation time will take longer (SS may then be
%		smaller than S).
%
% Usage:
%  [fS,fH]=fusionma(C,S)
%
% Input parameters:
%  C,S - Increasing column vectors. We want to know the conjugate at
%	S(j). C gives the slope of planar input points.
%
%
% Output parameter:
%  fS - Slopes supporting the epigraph. This is a subset of S. Column vector.
%  fH - Index at which the slope support the epigraph. Column vector.
%
%
% Being called by:
%  fusion - Choose either of fusionca.m or fusionma.m
%
% Example:
%  X=[-5:0.5:5]'
%  Y=X.^2
%  C=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
%  S=C-0.25*ones(size(C))
%  [fS,fH]=fusionma(C,S)

% The BUILT-IN sort function is faster than
%looking at each slope with a for loop
[Z I]=sort([C;S]);

I=find(I>size(C,1));
J=I(2:size(I,1))-I(1:size(I,1)-1);
K=J-ones(size(J));
L= cumsum(K)+ones(size(K));

H=[I(1) L'+I(1)-1]';
fS=S;fH=H;

% This (may be) obscure computation ONLY aims at using Matlab syntax
% to do exactly the same as fusionca.m whose programming sticks to the
% LLT as written in ["Faster than the Fast Legendre Transform, the
% Linear-time Legendre Transform", Y. Lucet]
