function display_row(U, s, theta, nbr_row, nbr_col, num_row)

[N,P] = size(U);

% la matrice de convolution
B = zeros(N);
B(1,1) = -4*s;
B(2,1) = s;
B(N,1) = s;
B(1,2) = s;
B(1,N) = s;

% precompute the matrix
Bf = real( fft2(B) );
Uf = fft2(U);
Cf = ( (1-theta)*Bf + ones(N) ) ./ ( ones(N) - theta*Bf  );

for i=1:nbr_col
	Uf = Uf .* Cf;
	U = real( ifft2( Uf ) );
	Umax = max(max(U));
	U = U./Umax*255;

	SUBPLOT(nbr_row,nbr_col, i + nbr_col*(num_row-1));
	str = sprintf('t=%0.5g',i*s/N^2);

	image(U);
	axis image;
	% axis off;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
	title(str);
    if i==1
        ylabel( sprintf('\\theta=%.1f',theta) );
    end
end;
