function [C,d,z0,xlim,ylim,domain] = load_correlation(name,normalize)

% load_correlation - load correlation function from a list of typical examples
%
%   [C,d,z0,xlim,ylim] = load_correlation(name,normalize);
%
%   C(x,y,x1,y1) is a symbolic symmetric correlation function.
%   d is the dimensionality of the problem.
%   set normalize=1 to impose that the covariance is 1 on the diagonal.
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

domain = [];

% correlation function, should be symmetric
syms x y x1 y1;
z0 = [0;0]; % center point
xlim = [-1;1]; ylim = [-1;1];
switch name
    case 'gaussian1d'
        %% convolution against a Gaussian %%
        sigma = .07; % std
        C(x,y,x1,y1) = exp(- ( (x-x1)^2 )/(2*sigma^2) );
        d = 1;
    case 'laplace1d'
        %% Laplace transform in 1-D %%
        s = 4; % offset on covariance
        C(x,y,x1,y1) = 1/(x+x1);
        d = 1;
        z0 = [3;0]; % center point
        xlim = [0;10];
    case 'neuro-like-1d'
        % measuring along a line x=0
        % phi(x)= [  1/( x^2 + t^2 ) ]_{t in R}
        C(x,y,x1,y1) = 1/(x*x1*(x+x1));
        d = 1;
        z0 = [2;0]; % center point
        xlim = [.1;4];
        
    case 'neuro-like-1d-2'
        % measuring along two lines at x=0 and x=1
        % phi(x)= [  1/( x^2 + t^2 ),  1/( (1-x)^2 + s^2 ) ]_{s,t in R}
        C(x,y,x1,y1) = 1/(x*x1*(x+x1)) + 1/((1-x)*(1-x1)*(2-x-x1));
        d = 1;
        % center point,
        z0 = [.29;0];  % near 0.5 gives degeneracy
        xlim = [.1;.9]; 
                  
    case 'neuro-like-2d'
        % measuring along a line x=0
        % phi(x)= [  1/( x^2 + (y-t)^2 ) ]_{t in R}
        C(x,y,x1,y1) = (x+x1) / ( (x*x1) * ( (y-y1)^2 + (x+x1)^2 ) );
        z0 = [.5;0]; % center point
        xlim = [.1;1];
        ylim = [-4;4];
        d = 2;
        
    case 'neuro-like-box'
        % measuring along 4 lines line x=0, x=1, y=0, y1
        C0(x,y,x1,y1) = (x+x1) / ( (x*x1) * ( (y-y1)^2 + (x+x1)^2 ) );
        C(x,y,x1,y1) = C0(x,y,x1,y1) + C0(1-x,y,1-x1,y1) + C0(y,x,y1,x1) + C0(1-y,x,1-y1,x1);
        z0 = [.35;.2]; % center point
        z0 = [.4;.4]; % center point
        rho = .01;
        xlim = [rho;1-rho];
        ylim = xlim;
        d = 2;
        
    case 'neuro-like-disc'
        z0 = [.4; .3];
        xlim = [-1;1];
        ylim = xlim;        
        C0(x,y,x1,y1) = 1-(x^2+y^2)*(x1^2+y1^2);
        C1(x,y,x1,y1) = (1-x1^2-y1^2)*(1-x^2-y^2)*( (1-x*x1-y*y1)^2+(y*x1-x*y1)^2 );
        C(x,y,x1,y1) = 2*pi*C0(x,y,x1,y1)/C1(x,y,x1,y1);
        d=2;
        domain = @(x,y)x.^2+y.^2<=.99;
        
    case 'laplace1d-normalized'
        %% Laplace transform in 1-D, normalized in L^2 norm %%
        s = 4; % offset on covariance
        C(x,y,x1,y1) = 2*sqrt(x*x1)/(x+x1);
        d = 1;
        z0 = [3;0]; % center point
        xlim = [0.05;10];
        
    case 'isolaplace'
        %% 2D isotropic laplace-like phi(x) = [exp(-t*|x|^2)]_{t>0}%%
        C(x,y,x1,y1) = 1/( x^2+y^2+x1^2+y1^2 );        
        z0 = [0;4];
        xlim = z0(1) + [-1;1];
        ylim = z0(2) + [-1;1];
        d = 2;
        
    case 'gaussian2d'
        %% convolution against a Gaussian %%
        sigma = .5; % std
        C(x,y,x1,y1) = exp(- ( (x-x1)^2+(y-y1)^2 )/(2*sigma^2) );
        d = 2;     
        
    case 'gmixture'
        %% evaluation of x=mean and y=cov of 1D Gaussians 
        %    phi(x,y) : t -> 1/sqrt(y) * exp(-(x-t)^2/(2*y) ) 
        C(x,y,x1,y1) = 1/sqrt(y+y1) * exp(- ( (x-x1)^2 )/(2*( y+y1 )) );
        d = 2;
        z0 = [0;2^2];
        xlim = [-3;3];
        ylim = [.5;6].^2;  
        
    case 'gmixture2'
        % alternate parameterization
        C(x,y,x1,y1) = 1/sqrt(y^2+y1^2) * exp(- ( (x-x1)^2 )/(2*( y^2+y1^2 )) );
        d = 2;
        z0 = [0;2];
        xlim = [-3;3];
        ylim = [.5;6];
        
    case 'gmixture3'
        % alternate parameterization
        C(x,y,x1,y1) = 1/sqrt(y^3+y1^3) * exp(- ( (x-x1)^2 )/(2*( y^3+y1^3 )) );
        d = 2;
        z0 = [0;2^(2/3)];
        xlim = [-3;3];
        ylim = [.5;6].^(2/3);
        
    case 'gausslaplace'
        % Gaussian in X and Laplace in Y.
        s = 1; % width of the gaussian
        C(x,y,x1,y1) = exp(- (x-x1)^2 /(2*s^2) ) / (y+y1);
        d = 2;
        z0 = [0;1.5];
        xlim = [-3;3];
        ylim = [.1 5];
        
    case 'covindep'
        %% evaluation of covariance diag(x,y) of 2D Gaussians 
        %    phi(x,y) : (s,t) -> 1/sqrt(x*y) * exp(-s^2/(2*x)-t^2/(2*y) ) 
        C(x,y,x1,y1) = 1/sqrt( (x+x1)*(y+y1) );
        d = 2;
        z0 = [2;3];
        xlim = z0(1)+[-1;1];
        ylim = z0(2)+[-1;1];
    otherwise 
        error(['Unknown: ' name '.']);
end

if normalize
    % normalize the columns of Phi in L^2 norm, so covariance is 1 along
    % the diagonal
    C(x,y,x1,y1) = C(x,y,x1,y1)/sqrt( C(x,y,x,y) * C(x1,y1,x1,y1) );
end

end