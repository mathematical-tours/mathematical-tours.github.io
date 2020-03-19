% test for how convergence rate depends on anisotropy
alist = 1.001 + 100*linspace(0,1,4000).^4;

E = [];

niter = 1000;
gmode = 'search';

RateF = [];
RateE = [];
for it=1:length(alist)
    
    a = alist(it);
    x = .9; y = .3;
    
    F = []; E = [];
    for i=1:niter
        % min_t (x-t*gx)^2+a*(y-t*gy)^2
        %    (x-t*x)^2+a*(y-t*a*y)^2
        %    x^2*(1-t)^2+a*y^2*(1-a*t)^2
        %    x^2*(1-t)+y^2*a^2*(1-a*t)=0
        %   t*(x^2+a^3*y^2) = x^2+a^2*y^2
        %   t = (x^2+a^2*y^2)/(x^2+a^3*y^2)
        switch gmode
            case 'search'
                r = (x^2+a^2*y^2)/(x^2+a^3*y^2);
            case 'low'
                r = .1;
            case 'large'
                r = .52;
        end
        x = x - r*x;
        y = y - r*a*y;
        F(end+1) = ( x^2+a*y^2 )/2;
        E(end+1) = abs(x+1i*y);
    end
    m = max(find(E>1e-60));
    RateE(end+1) = exp( log(E(m)/E(1)) / m );
    m = max(find(F>1e-60));
    RateF(end+1) = exp( log(F(m)/F(1)) / m );
end

[~,i] = max(RateE); 
[~,j] = max(RateF); 
fprintf('Worse rate: iterate=%.3f, function=%.3f.\n', alist(i), alist(j));

clf;
subplot(2,1,1); hold on;
plot(alist, RateE,'LineWidth', 2);
title('Rate on the iterates');
axis tight; box on;
xlabel('Condition number'); ylabel('rate');
axis([min(alist) max(alist) 0 .65]);
subplot(2,1,2); hold on;
plot(alist, RateF,'LineWidth', 2);
% plot(alist, 1./sqrt(alist));
title('Rate on the function');
axis tight; box on;
xlabel('Condition number'); ylabel('rate');
axis([min(alist) max(alist) 0 .65]);
saveas(gcf, 'rates.png');


% plot(alist,log10(E));