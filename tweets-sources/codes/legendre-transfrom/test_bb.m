clc, home
echo on
%
% Graphical test of the bb.m function which computes the planar convex
% hull of a set of points.
%
% Test-script for the Linear-time Legendre Transform (Yves Lucet - 09/05/97)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off

global m0;comp=1;
clear Temps Fl;
  echo on
  %% press any key after a line is displayed to continue
  echo off


disp('    n  &  CPU bb & flops bb   ');

for p=[10:10:100],	%[100:200:5000]
  X=[-2:4/(p-1):2]';
  Y=p*randn(p,1);% Take a random vector

  F=flops;t=cputime;
   [l1 l2]=bb(X,Y);
  CPU=cputime-t;F=flops-F;

  Temps(comp)=CPU;
  Fl(comp)=F;
  comp=comp+1;

  figure(2);clf;
  plot(X,Y);hold on;plot(X,Y,'o');
  plot(l1,l2);plot(l1,l2,'+');
  title('Convex hull of random points');

  disp([sprintf('%5.0f  & %7.2f & %10.2f & %2.8f',p,CPU,F)]);
  pause
end;

hold off;figure(1);clf;
subplot(1,2,1);plot(Temps,'+');hold on;plot(Temps);title('BB CPU time');
subplot(1,2,2);plot(Fl,'+'); hold on;plot(Fl);title('BB Flops');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
