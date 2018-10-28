function vx = InterpVF(Xi,Yi,v,x)

vx = [];
for d=1:2
    vx(:,d) = interp2(Yi,Xi,v(:,:,d),x(:,2),x(:,1));
end
    
end