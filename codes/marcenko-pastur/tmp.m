
% In Random Matrix Theory, MP law gives the probability density function
% of singular values of large rectangular random matrices;
% when the dimensions of matrix tend to infinity.

% This contribution illustrates the PDF of matrix Y(N,N)=(T^-1)X*X^T, 
% where X is random matrix whose entries X_i,j are independent 
% and identically distributed random variables with zero mean
% and variance s^2. The program is applicable for both uniform and random
% distributions.

% Ref :
% Marchenko,V. A., Pastur, L. A. (1967) "Distribution of eigenvalues for some sets of
% random matrices", Mat. Sb. (N.S.), 72(114):4, 507536

% (c) Youssef KHMOU, Applied Mathematics ,30 January,2015.

clear;
close all;
N=400;
T=700;
% Ratio of matrix dimensions
c=N/T;

% Sample
x=randn(N,T); % Normal distribution
%x=rand(N,T);  % Uniform distribution

s=std(x(:));

% spectral matrix
r=x*x'/T;

%eigenvalues
l=eig(r);
% Probability Density Function 
% number of points for measurement.
n=50;
% Boundaries 
a=(s^2)*(1-sqrt(c))^2;
b=(s^2)*(1+sqrt(c))^2;

[f,lambda]=hist(l,linspace(a,b,n));
% Normalization
f=f/sum(f);
% Theoretical pdf
ft=@(lambda,a,b,c) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
F=ft(lambda,a,b,c);
% Processing numerical pdf
F=F/sum(F);
F(isnan(F))=0;

% Results
figure;
h=bar(lambda,f);
set(h,'FaceColor',[.75 .75 .8]);
set(h,'LineWidth',0.25);
xlabel('Eigenvalue \lambda');
ylabel(' Probability Density Function f(\lambda)');
title(' MarchenkoPastur distribution');
lmin=min(l);
lmax=max(l);
%axis([-1 2*lmax 0 max(f)+max(f)/4]);

hold on;
plot(lambda,F,'g','LineWidth',2);
hold off;
