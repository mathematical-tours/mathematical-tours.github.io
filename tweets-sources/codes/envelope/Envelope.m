
rep = '../results/envelope/';
[~,~] = mkdir(rep);

% generate a curve
n = 501;
t = (0:n-1)'/n;

% helpers
% compute normal
normalize = @(c)c./abs(c);
switch boundm
    case 'per'
        Normal = @(c)1i*normalize( c([2:n 1]) - c([n 1:n-1]) );
        dp = @(c)c([1:end 1]);
    otherwise
        dp = @(c)c;
        Normal = @(c)1i*normalize( c([2:n n]) - c([1 1:n-1]) );        
end


% parabola
c = (2*t-1) + 2i*(2*t-1).^3;
c = (2*t-1) + 2i*(2*t-1).^1;
boundm = 'nper';

% ellipse
U = 1;

Ulist = .2:.2:1.8;
for U = Ulist
    R = @(x)sign(x).*abs(x).^U;
    c = R(cos(2*pi*t))+2i*R(sin(2*pi*t));
    boundm = 'per';
    
    cn = Normal(c);
    h = ComputeEnvelope(c,cn, boundm);
    
    clf; hold on;
    sub = 5;
    PlotLines(c(1:sub:end),cn(1:sub:end), 'k');
    plot(dp(c), 'b', 'LineWidth', 2);
    plot(dp(h), 'r', 'LineWidth', 2);
    axis equal;
    u = max( abs([real(h); imag(h); real(c); imag(c)]) );
    u = 2;
    axis(1.05*[-u u -u u]); axis off;
    saveas(gcf, [rep 'ellipelike-' num2str(10*U) '.eps'], 'epsc');
end

%% 
% Interpolation by rotating normals

U = 1;
    R = @(x)sign(x).*abs(x).^U;
c = R(cos(2*pi*t))+2i*R(sin(2*pi*t));
cn = Normal(c);
boundm = 'per';

Theta = linspace(0,1,9)*pi/2;
drawmode = 'overlay';
drawmode = 'new';
clf; hold on;
for i=1:length(Theta)
    f = (i-1)/(length(Theta)-1);
    theta = Theta(i);
    cn1 = exp(theta*1i)*cn;
    h = ComputeEnvelope(c,cn1, boundm);    
    if strcmp(drawmode, 'new')
    clf; hold on;
    end
    sub = 5;
    if strcmp(drawmode, 'new')
    PlotLines(c(1:sub:end),cn1(1:sub:end), 'k');
    end
    plot(dp(c), 'b', 'LineWidth', 2);
    plot(dp(h), 'Color',[1-f 0 f], 'LineWidth', 2);
    axis equal;
    u = max( abs([real(h); imag(h); real(c); imag(c)]) );
    u = 2;
    axis(1.05*[-u u -u u]); axis off;
    drawnow;
    if strcmp(drawmode, 'new')
    saveas(gcf, [rep 'rotellipsse-' num2str(i) '.eps'], 'epsc');
    end
end

    saveas(gcf, [rep 'rotellipsse-all.eps'], 'epsc');
    
return;

% quarter circle
n = 501;
t = linspace(0,1,n);
c = t; 
cn = (1-t)*1i - t;
h = ComputeEnvelope(c,cn, 'nper');
clf; hold on;
plot(h, 'r', 'LineWidth', 2);
PlotLines(c(1:20:end),cn(1:20:end), 'b');
axis equal; 
axis([0 1 0 1]); axis off;

