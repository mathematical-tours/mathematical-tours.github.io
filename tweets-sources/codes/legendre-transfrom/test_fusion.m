clc, home
echo on
%%
%% Test file for fusionca.m and fusionma.m
%% Comparison of fusionma.m with fusionca.m
%%
%% Test-script for the Linear-time Legendre Transform (Yves Lucet - 09/05/97)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo off

% Beware that taking random numbers affects the computation time.
% Indeed if s_1<...<s_m<c_1<...<c_p, both functions are very fast
% (since they do not look at all slopes S). Unfortunately, this case gives
% poor numerical accuracy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Temps Tempsca Fl Flca;

comp=1;disp('    q  &  CPU    & CPUc    &  flops_ma   & flops_ca   ');


%for p=[500:250:5000],
for p=[10:10:100],
%  C=sort(normrnd(1,1,p,1));

  C=sort(randn(p,1));
  % the choice S=sort(normrnd(1,1,p,1)); would only give an EXPECTED
  % linear computation time, so we do not take random slopes.
  S=(C(1:size(C,1)-1)+C(2:size(C,1)))./2;

  Fma=flops;tma=cputime;
   [l1 l2]=fusionma(C,S);
  CPUma=cputime-tma;Fma=flops-Fma;

  Fca=flops;tca=cputime;
   [l1ca l2ca]=fusionca(C,S);
  CPUca=cputime-tca;Fca=flops-Fca;

  Tempsma(comp)=CPUma;Flma(comp)=Fma;
  Tempsca(comp)=CPUca;Flca(comp)=Fca;
  comp=comp+1;
  disp([sprintf('%5.0f  & %7.2f & %7.2f & %10.0f & %10.0f ',p,CPUma,CPUca,Fma,Fca)]);
end;

hold off;clf;
n=[10:10:100];%n=[500:250:5000];
subplot(2,2,1);plot(n,Tempsma,'+');hold on;
plot(n,Tempsma);title('CPU time for fusionma');
subplot(2,2,2);plot(n,Flma,'+'); hold on;
plot(n,Flma);title('Flops for fusionma');

subplot(2,2,3);plot(n,Tempsca,'+');hold on;
plot(n,Tempsca);title('CPU time for fusionca');
subplot(2,2,4);plot(n,Flca,'+'); hold on;
plot(n,Flca);title('Flops for fusionca');
clear n;

echo on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From the test above, we deduce that fusionca takes less flops than
% fusionma, netherveless fusionma is faster because it is a matlab
% built-in function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off

