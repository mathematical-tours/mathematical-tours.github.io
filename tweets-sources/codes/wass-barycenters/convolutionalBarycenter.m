function [barycenter,v] = convolutionalBarycenter(p,alpha,areaWeights,kernel,kernelTranspose,entropyLimit, options)

% convolutionalBarycenter - compute optimal transport barycenter with regularization
%
%   [barycenter,v] =
%   convolutionalBarycenter(p,alpha,areaWeights,kernel,kernelTranspose,entropyLimit, options);
%
%   p is a (N,P) matrix where each p(:,i) is an input density.
%   barycenter is N vector output (the barycenter)
%   u is an (N,P) set of dual vectors (actually exp of the dual vectors
%
%   kernel takes as input an (N,P) matrix and blur each enty
%   kernelTranspose is the adjoint filtering (=[] assume kernel is
%       symmetric)
%   entropyLimit is a scalar imposing a fixed entropy lower bound 
%       to avoid over-smoothing
%       (=[] for no entropy bound)
%   areaWeights takes into account non-uniform grids (=[] will takes ones(N,1)).
%
%   Copyright (c) 2014 Justin Solomon

if nargin<5 || isempty(kernelTranspose)
    kernelTranspose = kernel; % assume symmetry
end
if isempty(areaWeights)
    areaWeights = ones(size(p,1),1); 
end

options.null = 0;
niter = getoptions(options, 'niter', 1500);
tol = getoptions(options, 'tol', 1e-7);
verb = getoptions(options, 'verb', 1);
displayFunction = getoptions(options, 'disp', []);
disp_rate = getoptions(options, 'disp_rate', 10);
v = getoptions(options,'initial_v',ones(size(p)));
barycenter = getoptions(options,'initial_barycenter',ones(size(p,1),1));
unit_area_projection = getoptions(options,'unit_area_projection',0);

alpha = alpha / sum(alpha);

for j=1:niter
    oldBarycenter = barycenter;

    w = p./(kernelTranspose(bsxfun(@times,v,areaWeights)));
    
    if unit_area_projection == 1
         integrals = sum(bsxfun(@times,areaWeights,v.*kernel(bsxfun(@times,w,areaWeights))),1);
         w = bsxfun(@rdivide,w,integrals);
    end
    
    d = v.*kernel(bsxfun(@times,w,areaWeights));

    d(d<1e-300) = 1e-300;

    % Log-sum-exp to multiply a bunch of things
    barycenter = exp(sum(bsxfun(@times,alpha,log(d)),2));

    entropy = -sum(areaWeights.*(barycenter.*log(barycenter)));

    % on the first iteration it fails?
    if j>1 && nargin>5 && ~isempty(entropyLimit) && entropy>entropyLimit
        % project onto a lower entropy value
        fn = @(x) full(-sum(x*areaWeights.*((barycenter.^x).*log(barycenter))) - entropyLimit);
        options = optimset('Display','none','tolfun',1e-4,'tolx',1e-4);
        try
            a = fzero(fn,[0.5 3],options); % never seen an 'a' outside this range
            if verb==1
                fprintf('\ta = %g\n',a);
            end
        catch
            a = 1;
            fprintf('\tProjection failed.\n');
        end
        barycenter = barycenter.^a;
    end

    v = bsxfun(@times,v,barycenter)./d;
    
    if unit_area_projection == 1
         integrals = sum(bsxfun(@times,areaWeights,v.*kernel(bsxfun(@times,w,areaWeights))),1);
         v = bsxfun(@rdivide,v,integrals);
    end
        
    if strcmp(class(displayFunction),'function_handle') && mod(j,disp_rate)==1
        displayFunction(barycenter);
        if ndims(barycenter)==1
            box on; set(gca, 'XTick', []); set(gca, 'YTick', []);
        end
        drawnow;
    end

    % Gabriel: I remove the "sqrt(" in the following formula 
    change = sum(abs(oldBarycenter-barycenter).*areaWeights);
    area = sum(barycenter.*areaWeights);
    
    if verb==1
        fprintf('Iteration %d:  change = %g, area = %g\n',j,full(change),full(area));
    elseif verb==2
        progressbar(j,niter);
    end    

    if j>2 && change<tol %&& abs(area-1) < 1e-5
        if verb==2
            progressbar(niter,niter);
        end
        return;
    end
end
