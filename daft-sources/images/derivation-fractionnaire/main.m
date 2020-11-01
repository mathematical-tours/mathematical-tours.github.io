clear;

N	= 256;

% pré-calcule de la Gaussienne et de sa transformée
ramp	= [-N:N-1];
ramp	= pi * ramp/N;
f	= exp( -(ramp.^2)/(0.5) );
f	= f - mean(f);	% recentrage de la gaussienne pour assurer continuité

% calcule les différentes dérivées fractionaires
i = 1;
for n = 0 : 2/15 :2
    fn = der_frac(f,n);
	fn 	= fn / max(abs(fn));		% renormalise
	% dessine dans une boite différente 
	subplot(4,4,i);
	t = sprintf( '%0.2f', n );	
	plot( real(fn) ); 
	text( 20,1.0, sprintf( '%0.2f', n ) );
	set( gca, 'Xtick', [], 'Ytick', [] );	
	axis( [0 2*N -1.2 1.2] ); axis square;
	h = line( [N N], [-1.2 1.2] );
	set( h, 'LineStyle', ':' );
	i = i + 1;
end