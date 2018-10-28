%%
% draw triangles as spirals

if not(exist('test'))
    test = 0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));


% input points
clf; hold on;
Z = []; Zlist = {}; Zold = [];
while true
    axis equal; axis([0 1 0 1]); 
    [a,b,button] = ginput(1);    
    if button==3
        if length(Z)<3
            break;
        else
            Zlist{end+1}=Z;
            plot(Z([1:end 1]), 'k');
            Z = [];
        end
    else
        Z(end+1) = a+1i*b;
        [s,k] = min(abs(Z(end)-Zold));
        if s<.05
            Z(end) = Zold(k);
        end
        Zold(end+1) = Z(end); 
    	plot(Z, '.b');
        plot(Z(end), 'r.', 'MarkerSize', 20);
    end
end


col = 'r';
lw = 1;

q = 50;
delta_list = linspace(.5,0.01,q);
delta_list = rescale(linspace(1,0,q).^3,.01,.4);


for it=1:q
    t = (it-1)/(q-1);
    delta = delta_list(it);
    clf;
    DrawTriangle(Zlist, delta, [t 0 1-t], lw);
    drawnow;
    saveas(gcf, [rep znum2str( it,2 ) '.png']);
end

% AutoCrop(rep, '');