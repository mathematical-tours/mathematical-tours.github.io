%%
% Display roots of polynomials

intmode = 'linear';
intmode = 'circle';

if not(exist('test'))
    test=0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));


% click and play
for k=1:2
    clf; hold on; 
    if k>1
        plot(z0{1}, 'ko');
    end
    z0{k} = [];
    while true
        axis equal; axis([-1 1 -1 1]);
        box on; set(gca,'XTick',[], 'YTick',[]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 15);
        if (k==1 & button==3)            
            break;
        end
        z0{k}(end+1) = a+1i*b;
        if (k==2) && (length(z0{1})==length(z0{2}))
            break;
        end
    end
    z0{k} = z0{k}(:);
    % build matrix
    [Q,R] = qr(randn(length(z0{k})));
    M{k} = Q*diag(z0{k})*Q';
end



Z = @(t)eig(M{1}*t + M{2}*(1-t));

q = 50; % #anim
r = q*5;

for i=1:q
    clf; hold on;
    for j=1:r
        t = (j-1)/(r-1);
        plot(Z(t), '.', 'MarkerSize', 30, 'Color', .8 + .2*[t 0 1-t]);
    end
    t = (i-1)/(q-1);
    plot(Z(t), '.', 'MarkerSize', 30, 'Color', [t 0 1-t]);
    axis equal; axis([-1 1 -1 1]);
    box on; set(gca,'XTick',[], 'YTick',[]);
    drawnow;
    saveas(gcf, [rep 'anim-' znum2str(i,2) '.png']);
end

% AutoCrop(rep, 'anim');