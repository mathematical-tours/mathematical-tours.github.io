function y = somme_glissante(x,P,k)
N = length(x); y =zeros(N,N);
for( u=N:-1:1 ) for( v=N:-1:1 )
    if( u<N ) y(u,v)=y(u,v)+y(u+1,v); end;
    if( v<N ) y(u,v)=y(u,v)+y(u,v+1); end;
    if( u<N && v<N ) y(u,v)=y(u,v)-y(u+1,v+1); end;
    y(u,v)=y(u,v)+x(u,v)^k;
    if( u+P<=N ) y(u,v)=y(u,v)-x(u+P,v)^k; end;
    if( v+P<=N ) y(u,v)=y(u,v)-x(u,v+P)^k; end;
    if( u+P<=N && v+P<=N ) y(u,v)=y(u,v)+x(u+P,v+P)^k; end;    
end; end;