%%
% display of evolution of Fourier descriptors



names = {'cat' 'cat'};
RS = [0 1];
names = {'cat' 'elephant'};
RS = [0 0];

addpath('../toolbox/');
rep = MkResRep(['descriptors/' names{1} '-' names{2}]);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20, 'XTick', [], 'YTick', []);



n0 = 512; % points on the shapes
% resample
p = 1024/4;
curvabs = @(gamma)[0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
resample1 = @(gamma,d)interp1(d/d(end),gamma,(0:p-1)'/p, 'linear');
resample = @(gamma)resample1( [gamma;gamma(1)], curvabs( [gamma;gamma(1)] ) );

CL = {};
for s=1:2
    %
    I = load_image(names{s}, n0);
    I = rescale(sum(I,3));
    %
    t = linspace(0,1,n0);
    C = contour(t,t,I,[.5,.5]);
    m = C(2,1);  c = C(:,2:m+1); c = c(1,:)'+1i*c(2,:)';
    c = resample(c);
    if RS(s)==1
       c = (c-.5-.5i)*.5*exp(2i*pi/3) + .5 +.5i ;
    end
    if s==2 && RS(s)==0
        % find optimal realignement
        e = [];
        for d=1:n0
            e(d)=norm(CL{1}-circshift(c,d));            
        end
        [~,d] = min(d);
        c = circshift(c,d);
    end
    %
    CL{s} = c;
end

m = 20; % #of Fourier descriptors
ex = @(x)x(2:m+1);
normalize = @(x)x/sum(x);
FD = @(c1)normalize(ex(abs(fft(c1))));

vmax = max(max(FD(CL{1})), max(FD(CL{1})));

q = 50; 
for it=1:q    
    t = (it-1)/(q-1); 
    c1 = CL{1}*(1-t)+CL{2}*t;
    %
    fd = FD(c1); 
    clf;
    bar(fd, 'EdgeColor', [t 0 1-t], 'FaceColor', [t 0 1-t]);
    axis([.5 m+.5 0 vmax*1.05]);
    SetAR(1/2)
    % drawnow;
    saveas(gcf, [rep 'descr-' znum2str(it,2) '.png']);
    clf;
    plot(c1([1:end 1]), 'Color',  [t 0 1-t], 'LineWidth', 2);
    axis equal; axis([0 1 0 1]); axis off; axis ij;
    drawnow;
   % saveas(gcf, [rep 'shape-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, ['descr-'])
% AutoCrop(rep, ['shape-'])
