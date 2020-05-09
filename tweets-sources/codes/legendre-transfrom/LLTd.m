function [SS,Conj]=LLTd(X,Y,S)
% LLTd		Compute numerically the discrete Legendre transform
%		of a set of planar points (X(i),Y(i)) at
%		slopes S, i.e., Conj(j)=max_{i}[S(j)*X(i)-Y(i)].
%		It uses the Linear-time Legendre Transform algorithm
%		to speed up computations.
%
%		All slopes not needed to build the affine interpolation may be
%		removed by using fusionca.m instead of fusionma.m in
%		fusion.m but computation time will take longer (SS may
%		then be smaller than S).
%
% Usage:
%  [SS,Conj]=LLTd(X,Y,S)
%
% Input parameters:
%  X,Y - Column vectors. Usually Y(i)=f(X(i))
%  S - Column vector. We want to know the conjugate at S(j)
%
%
% Output parameter:
%  SS - Slopes needed to define the conjugate. Removing affine parts
%	(by using fusionca.m instead of fusionma.m in fusion.m) may
%	result in needing less slopes. This is a column vector.
%  Conj - numerical values of the conjugate at slopes SS. This is
%	a column vector.
%
% Functions called:
%  bb - Compute the planar convex hull.
%  fusion - sort two increasing sequences of slopes. This is the core
%  of the LLT algorithm.
%
% Being called by:
%  LFt - Compute the conjugate of a function.
%
% Example:
%  X=[-5:0.5:5]'
%  Y=X.^2
%  S=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
%  [SS Conj]=LLTd(X,Y,S)


[XX,YY]=bb(X,Y); % Compute the convex hull (input sensitive algorithm).

h=size(XX,1);
C=(YY(2:h)-YY(1:h-1))./(XX(2:h)-XX(1:h-1)); % Compute the slopes associated
					    % with the primal points

[SS,H]=fusion(C,S);
% Merge sequences C and S. Edit fusion.m to emphasize computation
% speed (use fusionma.m) or to emphasize compression of output data (use
% fusionca.m).
%% j=H(i,1) is the first indice of the point X(j) such that the line
%% with slope S(i) supports the epigraph of the function at X(j).

Conj=SS.*XX(H(:,1)) -YY(H(:,1));

