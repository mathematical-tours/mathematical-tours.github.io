%%
% Non linear diffusions of porous type


rep = '../results/jko-flow/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


n = 1024/2;
p = 4;
H = @(u)u.^p;

t = (0:n-1)'/n;
f0 = max(1-( (t-.5)/.15 ).^2,0);

a = [2:n n]; b = [1 1:n-1];
Delta = @(f)(2*f - f(a) - f(b) )/2;

tau = .2;
niter = 5000*10;
if p==1
    niter = 50000;
end
if p==1.5
    niter = 100000;
end
if p==2
    niter = 150000;
end
if p==4
    niter = 10*100000;
end
ndisp = max(1,ceil(niter/20));

clf; hold on;
f = f0;
F = [];
for i=1:niter
    F(:,end+1) = f;
    if mod(i,ndisp)==1
        c = (i-1)/(niter-1);
        plot(f, 'LineWidth', 2, 'Color', [1-c,0,c]);
        drawnow;
    end
    f = f - tau * Delta( H(f) );
    f = max(f,0);
end



% svg run
k = 0;
E = norm(f-f0);
nextE = 0;
kdisp = 20;
%
clf; hold on;
f = f0;
k = 0;
for i=1:niter    
    if norm(f-f0)>=nextE   
        k = k+1;     
        c = min( (k-1)/(kdisp-1), 1);
        plot(f, 'LineWidth', 2, 'Color', [1-c,0,c]);
        drawnow;
        % update runtime
        nextE = norm(f-f0) + E/kdisp;
    end
    f = f - tau * Delta( H(f) );
    f = max(f,0);
end
axis tight;
SetAR(1/2); box on;
set(gca, 'XTick', [], 'YTick', []);
saveas(gcf, [rep 'nlin-diffus-p' num2str(round(10*p)) '.eps'], 'epsc');

