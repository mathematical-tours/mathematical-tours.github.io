%%
% Test for various Minkowski sums.
% https://en.wikipedia.org/wiki/Brunn%E2%80%93Minkowski_theorem

names = {'bunny', ''};

names = {'cat', 'horse'};
names = {'elephant', 'bunny'};

names = {'france', 'triskel'};
names = {'', ''};

name = 'france';
name = 'elephant';

addpath('../toolbox/');
rep = MkResRep(name);

n = 256*4;
t = linspace(0,1,n);
[Y,X] = meshgrid(t,t);
fillin = @(p)inpolygon(Y(:),X(:),real(p),imag(p));
fillin = @(p)reshape(fillin(p),[n n]);

% resample a curve using a fixed number of points
curvabs = @(c)[0;cumsum( 1e-5 + abs(c(1:end-1)-c(2:end)) )];
resample1 = @(c,d,p)interp1(d/d(end),c,(0:p-1)'/p, 'linear');
resample = @(c,p)resample1( [c;c(1)], curvabs( [c;c(1)] ),p );
    
% draw a polygons
if strcmp(name, '')
    p = DrawPoly();
else
    f = load_image(name,n); f = rescale(sum(f,3))>.5;
    [C,h] = contour(f,[.5,.5]);
    m = C(2,1);  c = C(:,2:m+1); p = c(1,:)'+1i*c(2,:)';
    p = conj( resample(p,200) );
end
p = p(:);
p = ( p-mean(p) );
p = (1+1i)/2 + .5*p/max(abs(p));


convol = @(x,y)real(ifft2(fft2(x).*fft2(y)));
convols = @(x,y)fftshift(convol(fftshift(x),fftshift(y)));
mink = @(A,B)double(convols(A,B)>1e-6);


convol = @(x,r)real(ifft2(fft2(x).^r));
convols = @(x,r)fftshift(convol(fftshift(x),r));
thresh = @(x)x>max(x(:))*1e-6;
mink = @(A,r)double(thresh(convols(A,r)));

scale = @(p,sc)(p-(1+1i)/2)*sc+(1+1i)/2;

q = 50; 
rlist = 1:6; % linspace(1,4,q);
q = length(rlist);

Area = [];
for it=1:q
    r = rlist(it);
    t = (it-1)/(q-1); 
    % rescale the second shape
    A = fillin( scale(p,1/r) ); A(end/2,end/2) = 1;
    %
    H = 1-mink(A,r);   
    Area(it)=sum(H(:)==0)/(n*n);
    % color
    f = cat(3, (1-H)*t + H*1, (1-H)*0 + H*1, (1-H)*(1-t) + H*1);
    clf;
    imageplot(f); axis xy;
    drawnow;
    imwrite(f(end:-1:1,:,:), [rep 'anim-' znum2str(it,2) '.png']);
end



% AutoCrop(rep, 'anim');

