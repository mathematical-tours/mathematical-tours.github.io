function x = ProxLP(y,p,tau)

if length(y(:))>1
    for i=1:length(y(:))
        x(i)=ProxLP(y(i),p,tau);
    end
    x = reshape(x,size(y));
    return
end

x = linspace(0,abs(y)*1.1,200);
[~,i] = min( 1/2*(x-abs(y)).^2 + tau*abs(x).^p );
x = x*sign(x(i));

end