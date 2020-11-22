%%
% Solve TV denoising on the dual

n = 128*4;

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(name,it)saveas(gcf, [rep name '-' znum2str(it,3) '.png']);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

x0 = zeros(n,1); 
I = round([.1 .4 .7 .8]*n);
x0(I) = [.6 .4 -1.3 .3];
x0 = cumsum(x0);
x0 = rescale(x0,-.8,.8);
y = x0 + .1*randn(n,1);

clf; hold on;
plot(x0, 'b');
plot(y, 'r');

% gradient
nabla = diag(ones(n,1)) - diag(ones(n-1,1),-1); 
nabla(1,end) = -1;
% primal 
%   min_x 1/2*|x-y|^2 + lambda*|x|_1
% dual
%   min_{|v|_inf <= lambda} 1/2*|nabla^T v-y|^2

q = 70;
niter = 15000;
niter = 15000*2;

tau = 1.8/norm(nabla)^2;
lambda = .1;
lambda = .3;


lambda = 1.9; niter = 60000;
lambda = 4;


lambda = .2; niter = 4000;
tau = .2/norm(nabla)^2;


ndisp = round(1+(niter-1)*linspace(0,1,q).^4);
ndisp = unique(ndisp);
q = length(ndisp);

save_mode = 1;

v = zeros(n,1); x = y; 
jt = 1;
for it=1:niter
    if ndisp(jt)==it
        clf;
        if save_mode==0
            subplot(2,1,1);
        end
        hold on;
        plot(y, 'r', 'LineWidth', 2);
        plot(x, 'b', 'LineWidth', 2);
        axis([1 n -1 1]); box on;
        if save_mode==0
            subplot(2,1,2);
        else
            set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'XTick', [], 'YTick', []);
            mysaveas('primal', jt);
            clf; hold on;
        end
        I = find(abs(abs(v/lambda)-1)<1e-3);
        hold on;
        plot(v/lambda, 'color', [0 .5 0], 'LineWidth', 2);
        plot([1 n],-[1 1], 'k--', 'LineWidth', 2);
        plot([1 n],+[1 1], 'k--', 'LineWidth', 2);
        plot(I, v(I)/lambda, 'k.', 'MarkerSize', 20);
        axis([1 n -1.05 1.05]); box on;        
        if save_mode==1
            set(gca, 'PlotBoxAspectRatio', [1 1/3 1], 'XTick', [], 'YTick', []);
            mysaveas('dual', jt);
        end
        drawnow;
        jt = jt+1;
    end
    %
    v = v - tau * nabla*( nabla'*v - y );
    v = min(max(v,-lambda),lambda); % projection
    E(it) = 1/2*norm(nabla'*v - y)^2;
    % primal variable
    x = y - nabla'* v;
    F(it) = 1/2*norm(x-y)^2 + lambda*norm(nabla*x,1);
    %
%    mysaveas(it);
end
