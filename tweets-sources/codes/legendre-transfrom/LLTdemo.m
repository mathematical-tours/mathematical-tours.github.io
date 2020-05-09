echo off
%%
%% Demo-script for the Linear-time Legendre Transform (Yves Lucet - 09/05/97)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, home

echo on

%% This is a script file to demonstrate the use of the LLT package
%% which computes the discrete Legendre Transform in linear time.
%% It shows how to use LFt.m, LLTd.m, LLTd2D.m, infconv.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%First Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test procedure for Figure 6. The conjugate u^* of u cannot be
% explicitly written with standard functions. So one has to use
% numerical computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo off
clear all;
n=1000;m=1000;
[S Cu]=LFt('LLTdemo_exa3',1e-6,10,n,-10,10,m);

hold off;clf;plot(S,Cu,'r');

echo on
pause  %% press any key to continue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%First Test bis %%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same computation can be done knowing only a sample of a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off
clear all;
n=1000;m=1000;

a=1e-6;b=10;c=-10;d=10;
X=[a:(b-a)/(n-1):b]'; % Generate the primal points
for i=1:size(X,1), Y(i)=feval('LLTdemo_exa3',X(i));end; % sample the function
Y=Y';
S=[c:(d-c)/(m-1):d]'; % Generate the slopes

[S Cu] = LLTd(X,Y,S);

hold on;plot(S,Cu,'b');

echo on
pause  %% press any key to continue



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Second Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test procedure for Figure 7: Two-dimensional computation of a
% nonsmooth conjugate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off

clear all;
a=1.3;comp=25;

X1=[-a:2*a/(comp-1):a]';
X2=X1;
S1=[LLTdemo_exa4_duxy(-a,0):(LLTdemo_exa4_duxy(a,0)-LLTdemo_exa4_duxy(-a,0))/(comp-1):LLTdemo_exa4_duxy(a,0)]';
S2=[LLTdemo_exa4_duxy(0,-a):(LLTdemo_exa4_duxy(0,a)-LLTdemo_exa4_duxy(0,-a))/(comp-1):LLTdemo_exa4_duxy(0,a)]';
n1=size(X1,1);
n2=size(X2,1);
m1=size(S1,1);
m2=size(S2,1);

for i1=1:n1,
  for i2=1:n2,
	    Y(i1,i2)=LLTdemo_exa4_uxy(X1(i1),X2(i2));
  end;
end;

Cu=LLTd2D(X1,X2,Y,S1,S2);

figure(1);hold off;clf;surf(S1,S2,Cu);view(-70,15);

echo on
pause  %% press any key to continue



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Third Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test procedure for Figure 8: Solution of a Hamilton--Jacobi equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo off

clear all;
global t;
n1=50;n2=50;m=50;nt=10;

X1=[-6:12/n1:6]';
i1=17;i2=33;

% We take more points in the neighborhood of the critical
% points: -2 and 2. In addition, X1 is increasing.
X1=[X1(1:i1-1);[X1(i1):(X1(i1+1)-X1(i1))/(n1-1):X1(i1+1)]';...
X1(i1+2:i2-1);[X1(i2):(X1(i2+1)-X1(i2))/(n1-1):X1(i2+1)]';X1(i2+2:n2)];

for i=1:size(X1,1),
  Y1(i,1)=feval('LLTdemo_exa5_f',X1(i));
end;

h=size(X1,1);
Pf=(Y1(2:h)-Y1(1:h-1))./(X1(2:h)-X1(1:h-1));

clear XX Ms;
T=[0.4:4/nt:4]';
compt=1;
for t=0.4:4/nt:4,
  disp(nt+1-compt);
  compt=compt+1;

  % Here the critical points are -t and t
  X2=[-t:2*t/(n2-1):t]';
  X2=[[X2(1):(X2(2)-X2(1))/(n2-1):X2(2)]';X2(3:n2-2);...
  [X2(n2-1):(X2(n2)-X2(n2-1))/(n2-1):X2(n2)]'];

  for i=1:size(X1,1),
	    Y2(i,1)=feval('LLTdemo_exa5_thetat',X2(i));
  end;

  h=size(X2,1);
  Pg=(Y2(2:h)-Y2(1:h-1))./(X2(2:h)-X2(1:h-1));

  % We build the slopes. S is strictly increasing
  S=sort([Pf;Pg]);
  S1(1)=S(1);
  for i=2:1:size(S,1),
	    if S(i)>S1(size(S1,1)),
	      S1(size(S1,1)+1,1)=S(i);
	    end;
  end;
  S=S1;

  Yr=infconvd(X1,Y1,X2,Y2,S);
  X=sort([X1;X2]);
  if t==T(1), Xt0=X;end;
  XX(:,compt)=X;
  Ms(:,compt)=Yr; % Store slice t
end;

% We deal with the first slice (t=0) to obtain the same number of points
% in Ms(:,1) and in Ms(:,i), i>1
for i=1:size(Xt0,1),
  Y1(i,1)=feval('LLTdemo_exa5_f',Xt0(i));
end;

XX(:,1)=Xt0;
Ms(:,1)=Y1;

figure(1);hold off;clf;
TT=ones(size(XX,1),1)*[0 T'];
surf(XX,TT,Ms);view(120,30);hold on;
title('F(x,t)');xlabel('x');ylabel('t');

echo on
pause  %% press any key to continue



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fourth Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test procedure for Figure 9: Deconvolution computation
% We compute Dec(x):=sup_y[g(y)-f(y-x)]=(g^*-f^*)^*(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off

clear all;
global a c;
n=1000;m=1000;Sb=100;

a=1;c=2;
[S Cf]=LFt('LLTdemo_exa6_Fdec',1,3,n,-Sb,Sb,m);
Xf=[1:2/n:3]';
for i=1:n+1,
  Yf(i)=LLTdemo_exa6_Fdec(Xf(i));
end;
Yf=Yf';

a=3;c=4;
[S Cg]=LFt('LLTdemo_exa6_Fdec',1,7,n,-Sb,Sb,m);
Xg=[1:6/n:7]';
for i=1:n+1,
  Yg(i)=LLTdemo_exa6_Fdec(Xg(i));
end;
Yg=Yg';

X=[0:4/n:4]';
C=Cg-Cf;
[X Dec]=LLTd(S,C,X);

figure(1);hold off;clf;plot(X,Dec);hold on;plot(Xf,Yf);
plot(Xg,Yg);text(6,-1,'g(x)');text(4.5,0,'h(x)');
text(3,0.25,'f(x)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
