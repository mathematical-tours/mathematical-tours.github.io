function [f,gx,gy] = Rosenb(x,y)

f = (x-1).^2 + 100*(y-x.^2).^2;
gx = 2*(x-1) + 400*x.*(x.^2-y);
gy = 200*(y-x.^2);


end