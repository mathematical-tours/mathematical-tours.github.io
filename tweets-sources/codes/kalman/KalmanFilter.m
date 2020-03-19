%%
% Test of Kalman filter


addpath('../toolbox/');
rep = MkResRep();

d = 2; % dimension of the model
m = 3; % dimension of measurements

% update matrix
theta = .03*pi;
F = .99*[cos(theta),sin(theta); -sin(theta),cos(theta)];
% measurements matrix
H = randn(m,d);

% covariance
Q = .06 * eye(d)*.02;
R = 0.1 * eye(m)*.02;

% initial state, unknown
x = [1;0];
hx = x;
P = zeros(d);

X = x;
hX = hx;

q = 70;
niter = 200;
disp_list = round(linspace(1,niter,q));
ndisp = 1;

for i=1:niter
    if i==disp_list(ndisp)
        clf; hold on;
        for j=1:i-1
            s = (j-1)/(niter-1);
            plot(X(1,j:j+1), X(2,j:j+1), 'color', (1-s)*[1 .3 .3] + s*[.3 0 0], 'LineWidth', 2);
            plot(hX(1,j:j+1), hX(2,j:j+1), 'color', (1-s)*[.3 .3 1] + s*[0 0 .3], 'LineWidth', 2);
        end
        axis equal; axis off; axis([-1 1 -1 1]*1.2);
        drawnow;
        saveas(gcf, [rep 'anim-' znum2str(ndisp,2) '.png']);
        ndisp = ndisp+1;
    end
    
    % noisy update, unknown
    x = F*x + sqrtm(Q)*randn(d,1);
    % noisy measurements
    z = H*x + sqrtm(R)*randn(m,1);
    % predict
    hx1 = F * hx;
    P1 = F*P*F' + Q;
    % update
    y = z - H* hx; %   residual
    S = H * P1 * H' + R;  %  residual cov
    K = P1 * H' * pinv(S);  % Optimal Kalman Gain
    hx = hx1 + K * y;
    P = (eye(d)- K * H )* P1;
    % save
    X(:,end+1) = x;
    hX(:,end+1) = hx;
end
