%%
% Test for pinsker


vmax = 1;
% Pinsker inequality.
for it=1:q
    s=Pint(it,:);
    A = reshape(abs( TV(s,S) ./ sqrt(2*KL( s, S )) ), [n n]);
    A(Sigma==1) = min(A(Sigma==1),vmax)/vmax; % for KL
    % turn into color image
    I = 1+floor(clamp(A)*(r-1));
    U = reshape(CM(I,:), [n n 3]);
    for i=1:3
        Us = U(:,:,i); Us(Sigma==0) = 1; U(:,:,i) = Us;
    end
    % display quantized colormap
    A(Sigma==0) = NaN;
    clf; hold on;
    imagesc(x,x,permute(U, [2 1 3]));
    contour(x,x,A',linspace(0,1,r), 'k');
    plot( s(1) + s(2)*exp(1i*pi/3), 'r.', 'MarkerSize', 25 );
    axis image; axis off;
    drawnow;
    mysaveas('pinsker', it);
end

return;


KL = @(x,y)sum(x.*log(x./y),2);
TV = @(x,y)sum(abs(x-y),2);

x = rand(100000,3);
y = rand(100000,3);
x = x./sum(x,2);
y = y./sum(y,2);

% TV<=sqrt(2*KL)
E = TV(x,y)./sqrt(2*KL(x,y));
hist(E,100);
[~,i] = max(E);
x(i,:)
y(i,:)

