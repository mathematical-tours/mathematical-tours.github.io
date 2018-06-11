function C = randCost( n,m,c )
% Matrix of random values in {1,...,c*min(n,m)}
%

    if nargin < 3
        c = 1;
    end

    mx = min(n,m);

    C = int32(randi(c*mx,[n,m]));


end

