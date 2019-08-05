%%
% Illustrate product and convolution of gaussian functions.


addpath('../toolbox/');
rep = MkResRep();


% generate a covariance
R = @(t)[cos(t) sin(t); -sin(t) cos(t)];
C = @(t,s)R(t)*diag([1 s])*R(-t);

eta = .3^2;
c0 = C(pi/3,eta);


daw_cov(c0, [1 0 0]);
saveas(gcf, [rep 'c-fixed.png']);

q = 50; 
for it=1:q
    
    s = (it-1)/(q-1);
    c1 = C(pi*s,eta);
    
        
        
    daw_cov(c1, [0 0 1]);
    saveas(gcf, [rep 'c-input-' znum2str(it,2) '.png']);
    
    
    c = (c0+c1)/2;
    daw_cov(c, [0 1 0]);
    saveas(gcf, [rep 'c-convol-' znum2str(it,2) '.png']);


    c = inv( (inv(c0)+inv(c1))/2 );
    daw_cov(c, [1 .8 0]);
    saveas(gcf, [rep 'c-prod-' znum2str(it,2) '.png']);
    drawnow;

end

