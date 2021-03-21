%%%
% plot zeta function on the critical line

addpath('../toolbox/');
rep = MkResRep();
mysaveas = @(it)saveas(gcf, [rep 'anim-' znum2str(it,3) '.png']);


rmax = 100;
n = 3000;
z = 1/2 + 1i*linspace(0,rmax,n);
Ze = zeta(z);

q = 120;
clf; hold on;
for it=1:q
    s1 = (it-1)/q;
    s2 = it/q;
    m1 = floor(s1*(n-1))+1;
    m2 = floor(s2*(n-1))+1;
    plot( Ze(m1:m2), 'color', [s1 0 1-s1],'LineWidth',3 );
    axis equal; 
    axis([min(real(Ze)) max(real(Ze)) min(imag(Ze)) max(imag(Ze))] + [-1 1 -1 1]*.15);
    set(gca, 'XTick', [], 'YTick', []); box on;
    drawnow;
    mysaveas(it);
end