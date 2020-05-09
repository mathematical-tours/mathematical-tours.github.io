function [SS,Conj]=LLTdirect(X,Y,S)
% LLTdirect	Compute numerically the discrete Legendre transform
%		of a set of planar points (X(i),Y(i)) at
%		slopes S, i.e., Conj(j)=max_{i}[S(j)*X(i)-Y(i)].
%		It uses straight computation for a quadratic-time
%		algorithm.
%
%
% Usage:
%  [SS,Conj]=LLTdirect(X,Y,S)
%
% Input parameters:
%  X,Y - Column vectors. Usually Y(i)=f(X(i))
%  S - Column vector. We want to know the conjugate at S(j)
%
%
% Output parameter:
%  SS - This column vector equals S. It is only for coherency with the
%	LLTd.m function.
%  Conj - numerical values of the conjugate at slopes SS. This is
%	a column vector.
%
%
% Example:
%  X=[-5:0.5:5]'
%  Y=X.^2
%  S=(Y(2:size(Y,1))-Y(1:size(Y,1)-1))./(X(2:size(X,1))-X(1:size(X,1)-1))
%  [SS Conj]=LLTdirect(X,Y,S)


SS=S; % For compatibility with the LLTd function only.
Conj = max(X * S' - Y* ones(size(S))' )';

