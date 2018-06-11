function [x,y] = PickDistributions()

col = {'r' 'g' 'b'};
q = 15; % #points sampled at each click
s = .025; % gaussian std
it = 0;
x = []; y = [];
clf; hold on;
while true
    axis equal; axis([0 1 0 1]);
    [a,b,B] = ginput(1);
    if B==1
        c = 1; yv = 1; 
    elseif B==3
        c = 2; yv =-1; 
    else
        break;
    end
    %
    z = a+1i*b + s*(randn(q,1)+1i*randn(q,1));
    plot(z, '.', 'MarkerSize', 15, 'color', col{c});
    x = [x; real(z), imag(z)];
    y = [y; ones(q,1)*yv];
end


end