function [lambda,eta] = ContractionRatio(K)

n = size(K,1);

if 1
eta = -Inf;
for k=1:n
    for l=1:n
        for i=1:n
        for j=1:n
            eta = max(eta,  ...
                K(i,k)*K(j,l)/( K(j,k)*K(i,l) ) ...
            );
        end
        end
    end
end

lambda = ( sqrt(eta)-1 ) / ( sqrt(eta)+1 );

else

% linear rate: theory
Qmax = zeros(n)-Inf;
Qmin = zeros(n)+Inf;
for k=1:n
    for l=1:n
        for i=1:n
            Qmax(k,l) = max(Qmax(k,l), K(i,k)/K(i,l));
            Qmin(k,l) = min(Qmin(k,l), K(i,k)/K(i,l));
        end
    end
end
eta = max(max(Qmax/Qmin));
lambda = ( sqrt(eta)-1 ) / ( sqrt(eta)+1 );

end

end