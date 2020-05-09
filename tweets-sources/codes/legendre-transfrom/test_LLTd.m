clc, home
echo on
%%
%% Test file for LFt.m, LLTd.m, and LLTd2D.
%% We show that both functions take linear time
%% for computing one-variable (see test 1) or two-variable
%% conjugate (see test 2).
%% Test 3 compars LLTd.m (our improve algorithm) with
%% LLTdirect.m (straight computation)
%%
%% Test-script for the Linear-time Legendre Transform (Yves Lucet - 09/05/97)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%First Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test procedure for Figure 4: Linear computation time of the LLT when
% the function u is univariate and n,m are in  [25,300]
% To generate Figure 4 in the article take m=[250:250:3000] and
% n=[250:250:3000] but this can take some time on some computers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo off
clear all;
disp('m');
comn=1;comm=1;

for m=[25:25:300],
  for n=[25:25:300],
	    t1=cputime;
	    LFt('LLTdemo_exa1',-10,10,n,-5,5,m);
	    CPUllt=cputime-t1;
	    T(comn,comm)=CPUllt;
	    comn=comn+1;
  end;
  disp((300-m)/25);
  comn=1;
  comm=comm+1;
end;

hold off;figure(1);clf;surf(T);view(60,30);
hold on;xlabel('n');ylabel('m');zlabel('CPU time in s.');

echo on
pause  %% press any key to continue


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Second Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test procedure for Figure 5: Linear computation time of the LLT when
% the function is bivariate and n_i=m_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo off
clear all;
a=10;
compteur=1;

echo on
n=500:250:3000; %use n=500:250:10000; to obtain Figure 5 in the article
echo off

for i=1:size(n,2),
  comp=round(sqrt(n(i)));
  disp(size(n,2)-i);
  X1=[-a:2*a/(comp-1):a]';
  X2=X1;
  S1=X1;S2=S1; % since uqua^*=uqua we can take S_i=X_i
  n1=size(X1,1);
  n2=size(X2,1);
  m1=size(S1,1);
  m2=size(S2,1);
  for i1=1:n1,
	      for i2=1:n2,
		Y(i1,i2)=LLTdemo_exa2([X1(i1);X2(i2)]);
	      end;
  end;

  t1=cputime;
  LLTd2D(X1,X2,Y,S1,S2);
  CPUllt=cputime-t1;
  T(compteur)=CPUllt;
  compteur=compteur+1;
end;

figure(1);hold off;clf;plot(n',T);hold on;plot(n',T,'+');
xlabel('n');ylabel('CPU in s.');

echo on
pause  %% press any key to continue


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Third Test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We compar straight computation as done with LLTdirect.m versus the
% Linear-time algorithm as done with LLTd.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo off
global m0 ctps cfps;

ctps=0;cfps=0;erreur=0;comn=1;comm=1;
clear Temps Fl Tempsc Flc;

disp('   m   & CPU LLT &CPUdirect&  flops LLT & flops direct');

for m=[20:20:60],
  for n=[20:20:100],
	X=[10/n:1/n:10]';
	Y=test_LLTd_u(X);
	S=-5+[10/m:1/m:10]';

	F=flops;t=cputime;
	 [S1 C1]=LLTd(X,Y,S); %% Computation of the conjugate
	CPU=cputime-t;F=flops-F;

	Fc=flops;tc=cputime;
	 [S2 C2]=LLTdirect(X,Y,S);%% Computation of the conjugate
	CPUc=cputime-tc;Fc=flops-Fc;

	Temps(comn,comm)=CPU;Fl(comn,comm)=F;
	Tempsc(comn,comm)=CPUc;Flc(comn,comm)=Fc;
	comn=comn+1;

	Er(comn,comm)=max(C1-C2);
  end;
  disp([sprintf('%5.0f  & %7.2f & %7.2f & %10.0f & %10.0f ',m,sum(Temps(:,comm)),sum(Tempsc(:,comm)),sum(Fl(:,comm)),sum(Flc(:,comm)))]);
  comn=1;comm=comm+1;
end;

disp([sprintf('Maximal error between both algorithms: %2.8g',max(max(Er)))]);

hold off;figure(1);clf;
subplot(2,2,1);surf(Temps);view(30,30);title('LLT CPU time');
subplot(2,2,2);surf(Fl);view(-20,30);title('LLT Flops');

subplot(2,2,3);surf(Tempsc);view(30,30);title('Straight computation CPU time');
subplot(2,2,4);surf(Flc);view(-20,30);title('Straight computation Flops');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

