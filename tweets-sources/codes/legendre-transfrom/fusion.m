function [fS,fH]=fusion(C,S)
% fusion	merge two increasing sequences.
%		fusion(C,S) gives the resulting slopes fS and
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
%  [fS,fH]=fusion(C,S)
%
% Input parameters:
%  C,S - Increasing column vectors. We want to know the conjugate at
%	S(j). C gives the slope of planar input points.
%
%
% Output parameter:
%  fS - Slopes supporting the epigraph. This is a subset of S.
%  fH - Index at which the slope support the epigraph.
%
%
% Functions called (choose one from the two):
%  fusionca - sorting algorithm programmed in matlab. It gives a
%	linear-time algorithm and gives only the minimum number of
%	points to built an affine interpolation of the conjugate. This
%	is the function describe in the article [Faster than the Fast
%	Legendre Transform, the Linear-time Legendre Transform, Y. Lucet].
%  fusionma - use matlab sort function (which is built-in) to speed up
%	computations but do not remove points belonging to affine
%	parts of the conjugate.
%
% Being called by:
%  LLTd - Compute the discrete Legendre transform.
%
% Example:
%  X=[-5:0.5:5]'
%  Y=X.^2
%  C=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
%  S=C-0.25*ones(size(C))
%  [fS,fH]=fusion(C,S)

[fS,fH]=fusionma(C,S);	%use fusionca for removing points belonging to
			% affine parts (output sensitive algorithm).

