

addpath('../toolbox/');
rep = MkResRep();

n = 1024;
t = linspace(0,1,n)';
normalize = @(x)x/sum(x);
G = @(m,s)exp(-(t-m).^2./(2*s.^2));

smin = .02;
smax = .07;
k = 4;
for l=1:2
    M{l} = rand(1,k);
    S{l} = smin+(smax-smin)*rand(1,k);
    A{l} = normalize( .5 + rand(1,k) );
end

% click and play
for l=1:2
    clf; hold on;
    axis([0 1 smin smax]);
    if l==2
        plot(M{1},S{1}, 'x', 'MarkerSize', 15); 
    end
    M{l} = []; S{l} = [];
    for r=1:k    
        [M{l}(end+1),S{l}(end+1),button] = ginput(1);
        plot(M{l}(end),S{l}(end), '.', 'MarkerSize', 15);    
    end
    if l==1
        A{l} = linspace(1,.2,k)';        
    else
        A{l} = A{1}(:) + .2 * linspace(.2,1,k)';
    end
    A{l} = A{l}(:)';
    M{l} = M{l}(:)';
    S{l} = S{l}(:)';
end

% Fused display
clf; hold on;
for l=1:2
    u = l-1;
    scatter( M{l},S{l}, (A{l}-.15)*800, [u 0 1-u], 'filled' );
end
axis([0 1 smin smax]);
box on;
set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'XTick', [], 'YTick', []);
saveas(gcf, [rep 'fused.png'], 'png');

    
q = 2;
for it=1:q
    u = (it-1)/(q-1);
    
    % interpolate
    m = (1-u)*M{1} + u * M{2};
    s = (1-u)*S{1} + u * S{2};
    a = (1-u)*A{1} + u * A{2};
    Z = n*normalize( sum(G(m,s) .* a,2) );
    
    % display spikes
    clf; 
    scatter( m,s, (a-.15)*800, [u 0 1-u], 'filled' );
    axis([0 1 smin smax]);
    box on;
    set(gca, 'PlotBoxAspectRatio', [1 2/3 1], 'XTick', [], 'YTick', []);
%    drawnow;
    saveas(gcf, [rep 'spikes-' znum2str(it,2) '.png'], 'png');
    
    % display curves
    clf;
%    plot(t,Z, 'LineWidth', 2, 'Color', [u 0 1-u]);
    h = area(t, Z, 'FaceColor', [u 0 1-u], 'EdgeColor', [u 0 1-u]);
    %h.FaceAlpha = 0.5;
    axis([0 1 0 5]); box on;
    set(gca, 'PlotBoxAspectRatio', [1 1/2 1], 'XTick', [], 'YTick', []);
    drawnow;
    saveas(gcf, [rep 'function-' znum2str(it,2) '.png'], 'png');
    
end
