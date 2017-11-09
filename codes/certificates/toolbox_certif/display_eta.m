function display_eta(Eta,d,Z,z0,x0,xlim,ylim, options)

% display_eta - display 1D and 2D certificate
%
%   display_eta(Eta,d,Z,z0,x0,xlim,ylim);
%
%   Eta(x,y) is a symbolic function.
%   d is the dimensionality of the problem.
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

options.null = 0;
domain = getoptions(options, 'domain', []);


switch d
    case 1
        Q = 201;
        Xa = linspace(xlim(1),xlim(2),Q);
        % eta = double( EtaW(Xa,Xa*0) );
        eta = fast_formula(Eta,Xa,Xa*0);
        hold on;
        plot(Xa,eta, 'LineWidth', 2);
        plot([xlim(1) xlim(2)], [1 1], 'r--');
        stem(z0(1),1, 'r.', 'MarkerSize', 25);
        axis tight; box on;
    case 2
        Q = 201;
        Xa = linspace(xlim(1),xlim(2),Q); 
        Ya = linspace(ylim(1),ylim(2),Q);
        [Y,X] = meshgrid(Ya,Xa);
        % eta = double( EtaW(X(:),Y(:)) );        
        % eta = reshape(eta, size(X));
        eta = fast_formula(Eta,X,Y);
        eta = real( eta );
        if not(isempty(domain))
            A = domain(X,Y); 
            eta(A==0) = NaN;
        end
        hold on;
        imagesc(Ya,Xa,eta);
        contour(Ya,Xa,eta, 20, 'k');
        cm1 = parula(256); 
        t = linspace(0,1,256)';
        cm2 = (1-t)*cm1(end,:)+t*[1 0 0];
        cm = [cm1;cm2]; colormap(cm);
%        caxis([min(eta(:)) max(eta(:))]);
        caxis([min(eta(:)) 2-min(eta(:))]);
        % fit within the box xlim x ylim
        if not(isempty(Z))
            Zc = Z-repmat(mean(Z,1), [size(Z,1) 1]);
            r = 10;
            for k=1:size(Z,1)
                plot( z0(2)+[0 r*Zc(k,2)], z0(1)+[0 r*Zc(k,1)], 'r--' );
            end
        end
        if not(isempty(z0)) && isempty(x0)
            plot(z0(2),z0(1), '.r', 'MarkerSize', 25);
        end
        for i=1:size(x0,1)
            plot(x0(i,2),x0(i,1), '.r', 'MarkerSize', 25);
        end
        %plot(Zc(:,2)/r+z0(2), Zc(:,1)/r+z0(1), '.r', 'MarkerSize', 25);
        axis([min(ylim) max(ylim) min(xlim) max(xlim)]); 
        axis square; axis on; box on;
end

end