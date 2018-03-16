addpath('../../matlab/'); 



for n=[1000 5000 10000 20000]
C = int32( round(rand(n)*1e6) );
tic; [rho,varrho,u,v] = hungarianLSAP(C); toc
end

disp(['min cost = ',num2str(sum(u)+sum(v))]);
