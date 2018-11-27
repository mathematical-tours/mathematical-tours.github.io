%%
% Test for various Minkowski sums.
% https://en.wikipedia.org/wiki/Brunn%E2%80%93Minkowski_theorem

names = {'bunny', ''};

names = {'cat', 'horse'};
names = {'elephant', 'bunny'};

names = {'france', 'triskel'};
names = {'', ''};

if not(exist('test'))
    test = 0;
end
test = test+1;

addpath('../toolbox/');
rep = MkResRep(num2str(test));

n = 256*2;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
fillin = @(p)inpolygon(Y(:),X(:),real(p),imag(p));
fillin = @(p)reshape(fillin(p),[n n]);

% resample a curve using a fixed number of points
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );
    
% draw a polygons
p = {}; F = {};
for i=1:2
    if isempty(names{i})
        p{i} = DrawPoly();
    else
        f = load_image(names{i},n); f = rescale(sum(f,3))>.5;
        [C,h] = contour(f,[.5,.5]);
        m = C(2,1);  c = C(:,2:m+1); p{i} = c(1,:)'+1i*c(2,:)';
        p{i} = conj( resample(p{i},200) );
    end
    p{i} = p{i}(:);
    p{i} = ( p{i}-mean(p{i}) ); 
    p{i} = (1+1i)/2 + .5*p{i}/max(abs(p{i})); 
end


convol = @(x,y)real(ifft2(fft2(x).*fft2(y)));
convols = @(x,y)fftshift(convol(fftshift(x),fftshift(y)));
mink = @(A,B)double(convols(A,B)>1e-6);


q = 50; 
Area = [];
for it=1:q
    t = (it-1)/(q-1);    
    % rescale the second shape
    A = fillin((1-t)*(p{1}-(1+1i)/2)+(1+1i)/2); A(end/2,end/2) = 1;
    B = fillin(t*(p{2}-(1+1i)/2)+(1+1i)/2); B(end/2,end/2) = 1;
    %
    H = 1-mink(A,B);   
    Area(it)=sum(H(:)==0)/(n*n);
    % color
    f = cat(3, (1-H)*t + H*1, (1-H)*0 + H*1, (1-H)*(1-t) + H*1);
    clf;
    imageplot(f); axis xy;
    drawnow;
    imwrite(f(end:-1:1,:,:), [rep 'anim-' znum2str(it,2) '.png']);
end


for it=1:q
    t = (it-1)/(q-1);  
    clf; hold on;
    plot(1:q, sqrt(Area), 'k', 'LineWidth', 2); 
    plot(it, sqrt(Area(it)), '.', 'Color', [t 0 1-t], 'MarkerSize', 50);
    axis([1 q min(sqrt(Area))-.01 max(sqrt(Area))+.01]);
    set(gca, 'XTick', [], 'YTick', [], 'PlotBoxAspectRatio', [1 1/3 1]); box on;
    drawnow;
    saveas(gcf, [rep 'volume-' znum2str(it,2) '.png']);
end

% AutoCrop(rep, 'volume');

