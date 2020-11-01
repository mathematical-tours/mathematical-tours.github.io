% dessine des graphes cycliques.

n = 11;

S{1} = [-1,1];
S{2} = [-2,2];
S{3} = [-3,3];

k = 0;
for n=[11 12]
for s = 1:3
    
    k = k+1;
    SS = S{s};
    subplot(2,3,k);
    [A,xy] = gen_cyclic_graph(n,SS);
    gplot( A,xy, 'k.-' );
    axis([-1,1,-1,1]);
    axis square;
    axis off;
    str = [];
    for i=1:length(SS)
        str = [str, sprintf('%d', SS(i))];
        if i~= length(SS)
            str = [str, ','];
        end
    end
    title( sprintf('G=Z/%dZ, S=[%s]', n, str) );
    
end
end