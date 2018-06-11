%% 
% Test for edge collapse.

name = 'venus';
name = 'moomoo';
name = 'duck';

addpath('../toolbox/');
rep = MkResRep(name);


options.name = name;
[V,F] = read_off([name '.off']);
n = size(V,2);

if strcmp(name,'duck')
    V(3,:) = -V(3,:);
end

myaxis = @()axis([min(V(1,:)) max(V(1,:)) min(V(2,:)) max(V(2,:)) min(V(3,:)) max(V(3,:))]*1.3);

% Display full quality.
clf;
plot_mesh(V,F,options);
shading faceted;
camlight;
view(147,-30);

q = 50;

m = floor( (size(F,2) - 4)/2 );
it_disp = round(linspace(1,m,q));
itk = 1;

F1 = F;
V1 = V;
for i=1:m-1
    edges = compute_edges(F1);
    D = V(:,edges(1,:)) - V(:,edges(2,:));
    D = sum(D.^2,1);
    [tmp,k] = min(D);
    e = edges(:,k);
    V1(:,e(1)) = mean( V1(:,e),2 );
    V1(:,e(2)) = Inf;
    F1(F1==e(2)) = e(1);
    a = sum( diff(sort(F1))==0 );
    F1(:,a>0) = [];
    if size(F1,2)<=5
        clf; plot_mesh(V1,F1,options);
        shading faceted;
        myaxis(); view(2);
        view(147,-30); % duck
        camlight;
        drawnow;
        saveas(gcf, [rep name '-' znum2str(itk,2) '.png']);
        break;
    end
    if i==it_disp(itk) 
        % display
        clf; plot_mesh(V1,F1,options);
        shading faceted;
        myaxis(); view(2);
         view(147,-30); % duck
        camlight;
        drawnow;
        saveas(gcf, [rep name '-' znum2str(itk,2) '.png']);
        itk = itk+1;
    end
end

% AutoCrop(rep, [name '-']); 