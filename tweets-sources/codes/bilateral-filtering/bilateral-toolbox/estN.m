function Nest=estN(sigmar,T,eps)
% sigmar        : width of range Gaussian
% T             : dynamic range of image is [0,2T]
% eps           : kernel approximation accuracy 
% Nest          : approximation order
if  sigmar > 70    
    N=10;
elseif sigmar < 5
    N=800;
else 
    lam=(T/sigmar)^2;
    p = log(exp(1)*lam);
    q = -lam - log(eps);
    t = q*exp(-1)/lam;
    W = t - t^2 + 1.5*t^3 - (8/3)*t^4; 
    N = min(max(q/W,10),300);
    if sigmar < 30
        for iter = 1:5  
            N = N - (N*log(N)-p*N-q)/(log(N)+1-p);
        end
    end 
end
Nest = ceil(N);
end
