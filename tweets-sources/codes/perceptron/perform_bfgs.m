function [f, R, info] = perform_bfgs(Grad, f, options)

% perform_bfgs - wrapper to HANSO code
%
%   [f, R, info] = perform_bfgs(Grad, f, options);
%
%   Grad should return (value, gradient)
%   f is an initialization
%   options.niter is the number of iterations.
%   options.bfgs_memory is the memory for the hessian bookeeping.
%   R is filled using options.report which takes as input (f,val).
%
%   Copyright (c) 2011 Gabriel Peyre

n = length(f);
pars.nvar = n;
pars.fgname = @(f,pars)Grad(f);
options.x0 = f;
options.nvec = getoptions(options, 'bfgs_memory', 20); % BFGS memory
options.maxit = getoptions(options, 'niter', 1000);
options.prtlevel = 0;
options.normtol = eps;
options.tol = eps;
options.verb = 1;
% options.report = @(x,val)struct('E', val, 'timing', cputime()-t);
[f, R, energy, d, H, iter, info] = bfgs(pars,options);
% 
% if (info == 7)
%     disp(['Not satisfying Wolf conditions!'])
% end

function [x, R, f, d, H, iter, info, X, G, w, fevalrec, xrec, Hrec] = bfgs(pars, options)
%BFGS The BFGS quasi-Newton minimization algorithm, Version 2.0, 2010
%   Basic call:[x, R, f, d] = bfgs(pars) 
%   Full call: [x, R, f, d, H, iter, info, X, G, w, fevalrec, xrec, Hrec] = bfgs(pars,options)
%   Input parameters
%    pars is a required struct, with two required fields
%      pars.nvar: the number of variables
%      pars.fgname: string giving the name of function (in single quotes) 
%         that returns the function and its gradient at a given input x, 
%         with call   [f,g] = fgtest(x,pars)  if pars.fgname is 'fgtest'.
%         Any data required to compute the function and gradient may be
%         encoded in other fields of pars.
%    options is an optional struct, with no required fields
%      options.x0: each column is a starting vector of variables
%          (default: empty)
%      options.nstart: number of starting vectors, generated randomly
%          if needed to augment those specified in options.x0
%          (default: 10 if options.x0 is not specified)
%      options.maxit: max number of iterations
%          (default 1000) (applies to each starting vector)
%      options.nvec: 0 for full BFGS matrix update, otherwise specifies 
%           number of vectors to save and use in the limited memory updates
%          (default: 0 if pars.nvar <= 100, otherwise 10)
%      options.H0: 
%          for full BFGS: initial inverse Hessian approximation
%           (must be positive definite, but this is not checked)
%          for limited memory BFGS: same, but applied every iteration
%           (must be sparse in this case)
%          (default: identity matrix, sparse in limited memory case)
%      options.scale: 
%          for full BFGS: 1 to scale H0 at first iteration, 0 otherwise
%          for limited memory BFGS: 1 to scale H0 every time, 0 otherwise
%          (default: 1)
%      options.ngrad: number of gradients willing to save and use in
%           solving QP to check optimality tolerance on smallest vector in
%           their convex hull; see also next two options 
%          (default: min(100, 2*pars.nvar, pars.nvar + 10)
%          (1 is recommended if and only if f is known to be smooth)
%      options.normtol: termination tolerance on d: smallest vector in
%           convex hull of up to options.ngrad gradients
%          (default: 1e-6) 
%      options.evaldist: the gradients used in the termination test
%           qualify only if they are evaluated at points  approximately 
%           within distance options.evaldist of x
%          (default: 1e-4) 
%      options.fvalquit: quit if f drops below this value 
%          (default: -inf) 
%      options.xnormquit: quit if norm(x) exceeds this value
%          (default: inf)
%      options.cpumax: quit if cpu time in secs exceeds this 
%          (default: inf) (applies to total running time)
%      options.strongwolfe: 0 for weak Wolfe line search (default)
%                           1 for strong Wolfe line search
%          (strong Wolfe line search is not recommended for use with
%           BFGS; it is very complicated and bad if f is nonsmooth;
%           however, it can be useful to simulate an exact line search)
%      options.wolfe1: first Wolfe line search parameter 
%          (ensuring sufficient decrease in function value, default: 0)
%          (should be > 0 in theory, but 0 is fine in practice)
%      options.wolfe2: second Wolfe line search parameter 
%          (ensuring algebraic increase (weak) or absolute decrease (strong)
%           in projected gradient, default: 0.5)
%          (important in theory and practice that this is not 0 or 1, 
%           except that it can be set to 0 if an exact line search is to be
%           simulated, using options.strongwolfe = 1)
%      options.quitLSfail: 1 if quit when line search fails, 0 otherwise
%          (default: 1, except if options.strongwolfe = 1 and
%           options.wolfe2 = 0, simulating exact line search)
%          (0 is potentially useful if f is not numerically continuous)
%      options.prtlevel: one of 0 (no printing), 1 (minimal), 2 (verbose)
%          (default: 1)
%
%   Output parameters: 
%    all return information on the runs for each starting vector
%    x: the final iterates
%    f: the final function values 
%    d: the final smallest vectors in the convex hull of the saved gradients 
%     at termination (the final gradient if options.ngrad == 1)
%    H: final BFGS inverse Hessian approximation matrices
%     (full BFGS update only, symmetrized so they are exactly symmetric)
%     (nan if limited memory updates were used)
%    iter: number of iterations
%    info: reason for termination:
%     0: tolerance on smallest vector in convex hull of saved gradients met
%     1: max number of iterations reached
%     2: f reached target value
%     3: norm(x) exceeded limit
%     4: cpu time exceeded limit
%     5: f is inf or nan at initial point
%     6: direction not a descent direction due to rounding error
%     7: line search bracketed minimizer but Wolfe conditions not satisfied
%     8: line search did not bracket minimizer: f may be unbounded below
%    X: iterates where saved gradients were evaluated (see below) 
%    G: saved gradients used for computation of smallest vector in convex hull 
%      of gradients at points near final x
%    w: weights giving the smallest vector in the convex hull of the saved
%      gradients
%    fevalrec: record of all function values evaluated in all line searches, 
%      including the final accepted values (nans if options.strongwolfe = 1)
%    xrec: record of all x iterates 
%    Hrec: record of all H iterates
%     (not symmetrized, may not be symmetric because of rounding error)
%    Note: if there is more than one starting vector, then:
%      f, iter, info are vectors of length options.nstart
%      x, d are matrices of size pars.nvar by options.nstart
%      H, X, G, w, xrec, Hrec are cell arrays of length options.nstart, and 
%      fevalrec is a cell array of cell arrays
%    Thus, for example, d(:,i) = G{i}*w{i}, for i = 1,...,options.nstart
%
%   BFGS is normally used for optimizing smooth, not necessarily convex, 
%   functions, for which the convergence rate is generically superlinear.
%   But it also works very well for functions that are nonsmooth at their  
%   minimizers, typically with a linear convergence rate and a final 
%   inverse Hessian approximation that is very ill conditioned, as long 
%   as a weak Wolfe line search is used. This version of BFGS will work
%   well both for smooth and nonsmooth functions and has a stopping 
%   criterion that applies for both cases, described above.
%   See A.S. Lewis and M.L. Overton, Nonsmooth Optimization via BFGS, 2008.
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "bfgs".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter defaults
if nargin == 0
   error('bfgs: "pars" is a required input parameter')
end
if nargin == 1
   options = [];
end
options = setdefaults(pars, options);  % set most default options
options = setx0(pars, options); % augment options.x0 randomly
x0 = options.x0;
nstart = size(x0,2);
cpufinish = cputime + options.cpumax;
fvalquit = options.fvalquit;
xnormquit = options.xnormquit;
prtlevel = options.prtlevel;
% set other options
options = setdefaultsbfgs(pars, options); 
for run = 1:nstart
    if prtlevel > 0 & nstart > 1
        fprintf('bfgs: starting point %d\n', run)
    end
    options.cpumax = cpufinish - cputime; % time left
    if nargout > 10
        [x(:,run), R, f(run), d(:,run), HH, iter(run), info(run), X{run}, G{run}, w{run}, ...
           fevalrec{run}, xrec{run}, Hrec{run}] = bfgs1run(x0(:,run), pars, options);
    elseif nargout > 7 % avoid computing fevalrec, xrec, Hrec which are expensive as they grow inside the main loop
        [x(:,run), R, f(run), d(:,run), HH, iter(run), info(run), X{run}, G{run}, w{run}] = bfgs1run(x0(:,run), pars, options);
    else % avoid computing unnecessary cell arrays 
        [x(:,run), R, f(run), d(:,run), HH, iter(run), info(run)] = bfgs1run(x0(:,run), pars, options);
    end
    % make exactly symmetric (too expensive to do inside optimization loop}
    H{run} = (HH + HH')/2; 
    if cputime > cpufinish | f < fvalquit | norm(x) > xnormquit
        break
    end
end
if nstart == 1 % no point returning cell arrays of length 1
    H = H{1};
    if nargout > 10
        fevalrec = fevalrec{1};
        xrec = xrec{1};
        Hrec = Hrec{1};  % don't symmetrize
    end
    if nargout > 7
        X = X{1};
        G = G{1};
        w = w{1};
    end
end


function [x, R, f, d, H, iter, info, X, G, w, fevalrec, xrec, Hrec] = bfgs1run(x0, pars, options)
% Version 2.0, 2010
% make a single run of BFGS from one starting point
% intended to be called by bfgs.m
% outputs: 
%    x: final iterate
%    f: final function value
%    d: final smallest vector in convex hull of saved gradients
%    H: final inverse Hessian approximation
%    iter: number of iterations
%    info: reason for termination
%     0: tolerance on smallest vector in convex hull of saved gradients met
%     1: max number of iterations reached
%     2: f reached target value
%     3: norm(x) exceeded limit
%     4: cpu time exceeded limit
%     5: f or g is inf or nan at initial point
%     6: direction not a descent direction (because of rounding)
%     7: line search bracketed minimizer but Wolfe conditions not satisfied
%     8: line search did not bracket minimizer: f may be unbounded below 
%    X: iterates where saved gradients were evaluated
%    G: gradients evaluated at these points
%    w: weights defining convex combination d = G*w
%    fevalrec: record of all function evaluations in the line searches
%    xrec: record of x iterates
%    Hrec: record of H iterates
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "bfgs".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = pars.nvar;
fgname = pars.fgname;
normtol = options.normtol;
fvalquit = options.fvalquit;
xnormquit = options.xnormquit;
cpufinish = cputime + options.cpumax;
maxit = options.maxit;
nvec = options.nvec;
prtlevel = options.prtlevel;
strongwolfe = options.strongwolfe;
wolfe1 = options.wolfe1;
wolfe2 = options.wolfe2;
quitLSfail = options.quitLSfail;
ngrad = options.ngrad;
evaldist = options.evaldist;
verb = getoptions(options, 'verb', 1);
report = getoptions(options, 'report', @(x,v)v);
H0 = options.H0;
H = H0; % sparse for limited memory BFGS 
scale = options.scale;
x = x0;
if isstr(fgname)
    [f,g] = feval(fgname, x, pars);
else
    [f,g] = fgname(x, pars);
end
d = g;
G = g;
X = x; 
nG = 1;
w = 1;
dnorm = norm(g);
if nvec > 0 % limited memory BFGS
    S = [];
    Y = [];
end
iter = 0;
if nargout > 9
    % so outputs defined if quit immediately
    fevalrec{1} = nan; % cell array
    xrec = nan*ones(n,1); % not cell array
    Hrec{1} = nan; % cell array
end
if isnaninf(f) % better not to generate an error return
    if prtlevel > 0
        fprintf('bfgs: f is infinite or nan at initial iterate\n')
    end
    info = 5;
    return
elseif isnaninf(g)
    if prtlevel > 0
        fprintf('bfgs: gradient is infinite or nan at initial iterate\n')
    end
    info = 5;
    return
elseif dnorm < normtol
    if prtlevel > 0
        fprintf('bfgs: tolerance on gradient satisfied at initial iterate\n')
    end
    info = 0;
    return
elseif f < fvalquit
    if prtlevel > 0
        fprintf('bfgs: below target objective at initial iterate\n')
    end
    info = 2;
    return
elseif norm(x) > xnormquit
    if prtlevel > 0
        keyboard
        fprintf('bfgs: norm(x) exceeds specified limit at initial iterate\n')
    end
    info = 3;
    return
end
clear R;
for iter = 1:maxit
    if verb
        progressbar(iter,maxit);
    end
    R(iter) = report(x,f);
    if nvec == 0 % full BFGS
        p = -H*g;
    else % limited memory BFGS
        p = -hgprod(H, g, S, Y);  % not H0, as in previous version
    end
    gtp = g'*p;
    if gtp >= 0 | isnan(gtp) % in rare cases, H could contain nans
       if prtlevel > 0
          fprintf('bfgs: not descent direction, quit at iteration %d, f = %g, dnorm = %5.1e\n',...
              iter, f, dnorm)
       end
       info = 6;
       return
    end
    gprev = g;  % for BFGS update
    if strongwolfe 
        % strong Wolfe line search is not recommended
        % except to simulate an exact line search
        % function values are not returned, so set fevalrecline to nan
        fevalrecline = nan;
        [alpha, x, f, g, fail] = ...
            linesch_sw(x, f, g, p, pars, wolfe1, wolfe2, fvalquit, prtlevel);
        if wolfe2 == 0 % exact line search: increase alpha slightly to get 
                       % to other side of any discontinuity in nonsmooth case
            increase = 1e-8*(1 + alpha); % positive if alpha = 0
            x = x + increase*p;
            if prtlevel > 1
                fprintf(' exact line sch simulation: slightly increasing step from %g to %g\n', alpha, alpha + increase)
            end
            [f,g] = feval(pars.fgname, x, pars);
        end
    else % weak Wolfe line search is the default
        [alpha, x, f, g, fail, notused, notused2, fevalrecline] = ...
                linesch_ww(x, f, g, p, pars, wolfe1, wolfe2, fvalquit, prtlevel);
    end
    % for the optimality check:
    % discard the saved gradients iff the new point x is not sufficiently
    % close to the previous point and replace them by new gradient 
    if alpha*norm(p) > evaldist
        nG = 1;
        G = g;
        X = x;
    % otherwise add new gradient to set of saved gradients, 
    % discarding oldest if already have ngrad saved gradients
    elseif nG < ngrad
        nG = nG + 1;
        G =  [g G];
        X = [x X];
    else % nG = ngrad
        G = [g G(:,1:ngrad-1)];
        X = [x X(:,1:ngrad-1)];
    end
    % optimality check: compute smallest vector in convex hull of qualifying 
    % gradients: reduces to norm of latest gradient if ngrad == 1, and the
    % set must always have at least one gradient: could gain efficiency
    % here by updating previous QP solution
    if nG > 1
        [w,d] = qpspecial(G); % Anders Skajaa code for this special QP
    else
        w = 1;
        d = g; 
    end
    dnorm = norm(d);
    if nargout > 9
        xrec(:,iter) = x;
        fevalrec{iter} = fevalrecline; % function vals computed in line search
        Hrec{iter} = H;
    end
    if prtlevel > 1
        nfeval = length(fevalrecline);
        fprintf('bfgs: iter %d: nfevals = %d, step = %5.1e, f = %g, nG = %d, dnorm = %5.1e\n', ...
            iter, nfeval, alpha, f, nG, dnorm)
    end
    if f < fvalquit % this is checked inside the line search
        if prtlevel > 0
            fprintf('bfgs: reached target objective, quit at iteration %d \n', iter)
        end
        info = 2;
        return
    elseif norm(x) > xnormquit % this is not checked inside the line search
        if prtlevel > 0
            fprintf('bfgs: norm(x) exceeds specified limit, quit at iteration %d \n', iter)
        end
        info = 3;
        return
    end
    if fail == 1 % line search failed (Wolfe conditions not both satisfied)
        if ~quitLSfail
            if prtlevel > 1
                fprintf('bfgs: continue although line search failed\n')
            end
        else % quit since line search failed
            if prtlevel > 0
                fprintf('bfgs: quit at iteration %d, f = %g, dnorm = %5.1e\n', iter, f, dnorm)
            end
            info = 7;
            return
        end
    elseif fail == -1 % function apparently unbounded below
        if prtlevel > 0
           fprintf('bfgs: f may be unbounded below, quit at iteration %d, f = %g\n', iter, f)
        end
        info = 8;
        return
    end
    if dnorm <= normtol
        if prtlevel > 0 
            if nG == 1 
                fprintf('bfgs: gradient norm below tolerance, quit at iteration %d, f = %g\n', iter, f')
            else
                fprintf('bfgs: norm of smallest vector in convex hull of gradients below tolerance, quit at iteration %d, f = %g\n', iter, f')
            end
        end
        info = 0;
        return
    end
    if cputime > cpufinish
        if prtlevel > 0
            fprintf('bfgs: cpu time limit exceeded, quit at iteration %d\n', iter)
        end
        info = 4;
        return
    end
    s = alpha*p;
    y = g - gprev;
    sty = s'*y;    % successful line search ensures this is positive
    if nvec == 0   % perform rank two BFGS update to the inverse Hessian H
        if sty > 0 
            if iter == 1 & scale
                % for full BFGS, Nocedal and Wright recommend scaling I before 
                % the first update only
                H = (sty/(y'*y))*H; 
            end
            % for formula, see Nocedal and Wright's book
            rho = 1/sty;
            rhoHyst = rho*(H*y)*s';                                       % M = I - rho*s*y';
            H = H - rhoHyst' - rhoHyst + rho*s*(y'*rhoHyst) + rho*s*s';   % H = M*H*M' + rho*s*s';
        else % should not happen unless line search fails, and in that case should normally have quit
            if prtlevel > 1
                fprintf('bfgs: sty <= 0, skipping BFGS update at iteration %d \n', iter)
            end
        end
    else % save s and y vectors for limited memory update
        s = alpha*p;
        y = g - gprev;
        if iter <= nvec
            S = [S s];
            Y = [Y y];
        else % could be more efficient here by avoiding moving the columns
            S = [S(:,2:nvec) s];
            Y = [Y(:,2:nvec) y];
        end
        if scale 
            H = ((s'*y)/(y'*y))*H0;  % recommended by Nocedal-Wright
        end
    end
end % for loop
if prtlevel > 0
    fprintf('bfgs: %d iterations reached, f = %g, dnorm = %5.1e\n', maxit, f, dnorm)
end
info = 1; % quit since max iterations reached


function [xbundle, gbundle] = getbundle(x, g, samprad, N, pars);
%  get bundle of N-1 gradients at points near x, in addition to g,
%  which is gradient at x and goes in first column
%  intended to be called by gradsampfixed
%
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "gradsamp".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = length(x);
xbundle(:,1) = x;
gbundle(:,1) = g;
for k = 2:N  % note the 2
   xpert = x + samprad*(rand(m,1) - 0.5); % uniform distribution
   [f,grad] = feval(pars.fgname, xpert, pars);
   count = 0;
   while isnaninf(f) | isnaninf(grad)  % in particular, disallow infinite function values
       xpert = (x + xpert)/2;     % contract back until feasible
       [f,grad] = feval(pars.fgname, xpert, pars);
       count = count + 1;
       if count > 100 % should never happen, but just in case
           error('too many contractions needed to find finite f and grad values')
       end
   end; % discard function values
   xbundle(:,k) = xpert;
   gbundle(:,k) = grad;   
end


function options = setdefaults(pars,options)
%  call: options = setdefaults(pars,options)
%  check that fields of pars and options are set correctly and
%  set basic default values for options that are common to various 
%  optimization methods, including bfgs and gradsamp
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "bfgs".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    options = [];
end
if ~isfield(pars, 'nvar')
   error('setdefaults: input "pars" must have a field "nvar" (number of variables)')
elseif ~isposint(pars.nvar)
   error('setdefaults: input "pars.nvar" (number of variables) must be a positive integer')
end
if ~isfield(pars, 'fgname')
   error('setdefaults: input "pars" must have a field "fgname" (name of m-file computing function and gradient)')
end
if isfield(options, 'maxit')
    if ~isnonnegint(options.maxit)
        error('setdefaults: input "options.maxit" must be a nonnegative integer')
    end
else
    options.maxit = 1000;
end
if isfield(options, 'normtol')
    if ~isposreal(options.normtol)
        error('setdefaults: input "options.normtol" must be a positive real scalar')
    end
else
    options.normtol = 1.0e-6;
end
if isfield(options, 'fvalquit')
    if ~isreal(options.fvalquit)|~isscalar(options.fvalquit)
        error('setdefaults: input "options.fvalquit" must be a real scalar')
    end
else
    options.fvalquit = -inf;
end
if isfield(options, 'xnormquit')
    if ~isreal(options.xnormquit)|~isscalar(options.xnormquit)
        error('setdefaults: input "options.fvalquit" must be a real scalar')
    end
else
    options.xnormquit = inf;
end
if ~isfield(options, 'cpumax')
    options.cpumax = inf;
end
if ~isfield(options, 'prtlevel')
    options.prtlevel = 1;
end

function ipi = isposint(x)
% return true if x is a positive integer, false otherwise 
%
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isscalar(x) 
    ipi = 0;
else % following is OK since x is scalar
    ipi = isreal(x) & round(x) == x & x > 0;
end

function inni = isnonnegint(x)
% return true if x is a nonnegative integer, false otherwise
%
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isscalar(x) 
    inni = 0;
else % following is OK since x is scalar
    inni = (isreal(x) & round(x) == x & x >= 0);
end
 
function ipr = isposreal(x)
% return true if x is a positive real scalar, false otherwise 
%
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isscalar(x) 
    ipr = 0;
else % following is OK since x is scalar
    ipr = isreal(x) & x > 0;
end

function options = setx0(pars,options)
% set columns of options.x0 randomly if not provided by user
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "bfgs".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvar = pars.nvar;
if ~isfield(options, 'x0')
    options.x0 = [];
end
if isempty(options.x0)
    if isfield(options, 'nstart')
        if ~isposint(options.nstart)
            error('setx0: input "options.nstart" must be a positive integer when "options.x0" is not provided')
        else
            options.x0 = randn(nvar, options.nstart);
        end
    else
        options.x0 = randn(nvar, 10);
    end
else
    if size(options.x0,1) ~= nvar
        error('setx0: input "options.x0" must have "pars.nvar" rows')
    end
    if isfield(options, 'nstart')
        if ~isnonnegint(options.nstart)
            error('setx0: input "options.nstart" must be a nonnegative integer')
        elseif options.nstart < size(options.x0,2)
            error('setx0: "options.nstart" is less than number of columns of "options.x0"')
        else % augment vectors in options.x0 with randomly generated ones
            nrand = options.nstart - size(options.x0,2);
            options.x0 = [options.x0  randn(nvar, nrand)];
        end
    end % no else part, options.x0 is as provided
end

function options = setdefaultsbfgs(pars, options)
%  set defaults for BFGS (in addition to those already set by setdefaults)
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "bfgs".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line search options
if isfield(options, 'strongwolfe')
    if options.strongwolfe ~= 0 & options.strongwolfe ~= 1
        error('setdefaultsbfgs: input "options.strongwolfe" must be 0 or 1')
    end
else
    % strong Wolfe is very complicated and is bad for nonsmooth functions
    options.strongwolfe = 0;  
end
if isfield(options, 'wolfe1') % conventionally anything in (0,1), but 0 is OK while close to 1 is not
    if ~isreal(options.wolfe1) | options.wolfe1 < 0 | options.wolfe1 > 0.5
        error('setdefaultsbfgs: input "options.wolfe1" must be between 0 and 0.5')
    end
else
    options.wolfe1 = 1e-4; % conventionally this should be positive, and although
                           % zero is usually fine in practice, there are exceptions
end
if isfield(options, 'wolfe2') % conventionally should be > wolfe1, but both 0 OK for e.g. Shor
    if ~isreal(options.wolfe2) | options.wolfe2 < options.wolfe1  | options.wolfe2 >= 1
        error('setdefaultsbfgs: input "options.wolfe2" must be between max(0,options.wolfe1) and 1')
    end
else
    options.wolfe2 = 0.5;  % 0 and 1 are both bad choices
end
if options.strongwolfe
    if options.prtlevel > 0
        if options.wolfe2 > 0
            fprintf('setdefaultsbfgs: strong Wolfe line search selected, but weak Wolfe is usually preferable\n')
            fprintf('(especially if f is nonsmooth)\n')
        else
            fprintf('setdefaultsbfgs: simulating exact line search\n')
        end
    end
    if ~exist('linesch_sw')
        error('"linesch_sw" is not in path: it can be obtained from the NLCG distribution')
    end
else
    if ~exist('linesch_ww')
        error('"linesch_ww" is not in path: it is required for weak Wolfe line search')
    end
end
if isfield(options, 'quitLSfail')
    if options.quitLSfail ~= 0 & options.quitLSfail ~= 1
        error('setdefaultsbfgs: input "options.quitLSfail" must be 0 or 1')
    end
else
    if options.strongwolfe == 1 & options.wolfe2 == 0
        % simulated exact line search, so don't quit if it fails
        options.quitLSfail = 0;
    else  % quit if line search fails
        options.quitLSfail = 1;
    end
end
% other default options
n = pars.nvar;
if isfield(options, 'nvec')
    if ~isnonnegint(options.nvec)
        error('setdefaultsbfgs: input "options.nvec" must be a nonnegative integer')
    end
elseif n <= 100
    options.nvec = 0;  % full BFGS
else
    options.nvec = 10; % limited memory BFGS
end
if isfield(options,'H0')
    % H0 should be positive definite but too expensive to check
    if any(size(options.H0) ~= [n n])
        error('bfgs: input options.H0 must be matrix with order pars.nvar')
    end
    if options.nvec > 0 & ~issparse(options.H0)
        error('bfgs: input "options.H0" must be a sparse matrix when "options.nvec" is positive')
    end
else
    if options.nvec == 0
        options.H0 = eye(n); % identity for full BFGS
    else
        options.H0 = speye(n); % sparse identity for limited memory BFGS
    end
end
if isfield(options, 'scale')
    if options.scale ~= 0 & options.scale ~= 1
        error('setdefaultsbfgs: input "options.scale" must be 0 or 1')
    end
else
    options.scale = 1;
end
    
% note: if f is smooth, superlinear convergence will ensure that termination
% takes place before too many gradients are used in the QP optimality check
% so the optimality check will not be expensive in the smooth case
if isfield(options,'ngrad')
    if ~isnonnegint(options.ngrad)
        error('setdefaultsbfgs: input "options.ngrad" must be a nonnegative integer')
    end
else % note this could be more than options.nvec
     % rationale: it is only towards the end that we start accumulating
     % many gradients, and then they may be needed to veryify optimality
    options.ngrad = min([100, 2*pars.nvar, pars.nvar + 10]);
end
if isfield(options,'evaldist') 
    if ~isposreal(options.ngrad)
        error('setdefaultsbfgs: input "options.evaldist" must be a positive real scalar')
    end
else
    options.evaldist = 1e-4; 
end

function [alpha, xalpha, falpha, gradalpha, fail, beta, gradbeta, fevalrec] = ...
          linesch_ww(x0, f0, grad0, d, pars, c1, c2, fvalquit, prtlevel)
% LINESCH_WW Line search enforcing weak Wolfe conditions, suitable
%            for minimizing both smooth and nonsmooth functions
%          Version 2.0 for HANSO 2.0
% call:  [alpha, xalpha, falpha, gradalpha, fail, beta, gradbeta, fevalrec] = ...
%         linesch_ww_mod(x0, f0, grad0, d, pars, c1, c2, fvalquit, prtlevel);
%  Input
%   x0:      intial point
%   f0:      function value at x0
%   grad0:   gradient at x0
%   d:       search direction  
%   pars:    a structure that specifies the function name as well
%            anything else that the user needs to access in programming the
%            function and gradient values
%        pars.fgname:  name of function that returns function and gradient
%            it expects as input only x and pars, a parameter structure 
%            it is invoked by: [f,g] = feval(fgname, x, pars)
%   c1: Wolfe parameter for the sufficient decrease condition 
%          f(x0 + t d) ** < ** f0 + c1*t*grad0'*d     (DEFAULT 0)
%   c2: Wolfe parameter for the WEAK condition on directional derivative
%          (grad f)(x0 + t d)'*d ** > ** c2*grad0'*d  (DEFAULT 0.5)
%        where 0 <= c1 <= c2 <= 1.
%        For usual convergence theory for smooth functions, normally one
%        requires 0 < c1 < c2 < 1, but c1=0 is fine in practice.
%        May want c1 = c2 = 0 for some nonsmooth optimization 
%        algorithms such as Shor or bundle, but not BFGS.
%        Setting c2=0 may interfere with superlinear convergence of
%        BFGS in smooth case.
%   fvalquit: quit immediately if f drops below this value, regardless
%        of the Wolfe conditions (default -inf)
%   prtlevel: 0 for no printing, 1 minimal (default), 2 verbose 
%
%  Output:
%   alpha:   steplength satisfying weak Wolfe conditions if one was found,
%             otherwise left end point of interval bracketing such a point
%             (possibly 0)
%   xalpha:  x0 + alpha*d
%   falpha:  f(x0 + alpha d)
%   gradalpha:(grad f)(x0 + alpha d)  
%   fail:    0 if both Wolfe conditions satisfied, or falpha < fvalquit
%            1 if one or both Wolfe conditions not satisfied but an
%               interval was found bracketing a point where both satisfied
%           -1 if no such interval was found, function may be unbounded below
%   beta:    same as alpha if it satisfies weak Wolfe conditions,
%             otherwise right end point of interval bracketing such a point
%             (inf if no such finite interval found)
%   gradbeta: (grad f)(x0 + beta d) (this is important for bundle methods)
%             (vector of nans if beta is inf)
%             
%   fevalrec:  record of function evaluations

% The weak Wolfe line search is far less complicated that the standard 
% strong Wolfe line search that is discussed in many texts. It appears
% to have no disadvantages compared to strong Wolfe when used with
% Newton or BFGS methods on smooth functions, and it is essential for 
% the application of BFGS or bundle to nonsmooth functions as done in HANSO.
% However, it is NOT recommended for use with conjugate gradient methods,
% which require a strong Wolfe line search for convergence guarantees.
% Weak Wolfe requires two conditions to be satisfied: sufficient decrease
% in the objective, and sufficient increase in the directional derivative
% (not reduction in its absolute value, as required by strong Wolfe).
%
% There are some subtleties for nonsmooth functions.  In the typical case
% that the directional derivative changes sign somewhere along d, it is
% no problem to satisfy the 2nd condition, but descent may not be possible
% if the change of sign takes place even when the step is tiny. In this
% case it is important to return the gradient corresponding to the positive 
% directional derivative even though descent was not obtained. On the other 
% hand, for some nonsmooth functions the function decrease is steady
% along the line until at some point it jumps to infinity, because an
% implicit constraint is violated.  In this case, the first condition is
% satisfied but the second is not. All cases are covered by returning
% the end points of an interval [alpha, beta] and returning the function 
% value at alpha, but the gradients at both alpha and beta. 
%
% The assertion that [alpha,beta] brackets a point satisfying the
% weak Wolfe conditions depends on an assumption that the function 
% f(x + td) is a continuous and piecewise continuously differentiable 
% function of t, and that in the unlikely event that f is evaluated at
% a point of discontinuity of the derivative, g'*d, where g is the 
% computed gradient, is either the left or right derivative at the point
% of discontinuity, or something in between these two values.
%
% For functions that are known to be nonsmooth, setting the second Wolfe
% parameter to zero makes sense, especially for a bundle method, and for
% the Shor R-algorithm, for which it is essential.  However, it's not
% a good idea for BFGS, as for smooth functions this may prevent superlinear 
% convergence, and it can even make trouble for BFGS on, e.g., 
% f(x) = x_1^2 + eps |x_2|, when eps is small.
%
% Line search quits immediately if f drops below fvalquit.
%
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "bfgs".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6  % check if the optional Wolfe parameters were passed
    c1 = 0; % not conventional, but seems OK.  See note at top.
end
if nargin < 7
    c2 = 0.5; % see note at top
end
if nargin < 8
    fvalquit = -inf;
end
if nargin < 9
    prtlevel = 1;
end
if (c1 < 0 | c1 > c2 | c2 > 1) & prtlevel > 0 % allows c1 = 0, c2 = 0 and c2 = 1
   fprintf('linesch_ww_mod: Wolfe parameters do not satisfy 0 <= c1 <= c2 <= 1\n')
end
fgname = pars.fgname;
alpha = 0;  % lower bound on steplength conditions
xalpha = x0;
falpha = f0;
gradalpha = grad0; % need to pass grad0, not grad0'*d, in case line search fails
beta = inf;  % upper bound on steplength satisfying weak Wolfe conditions
gradbeta = nan*ones(size(x0));
g0 = grad0'*d; 
if g0 >= 0
    % error('linesch_ww_mod: g0 is nonnegative, indicating d not a descent direction')
    fprintf('linesch_ww_mod: WARNING, not a descent direction\n')
end
dnorm = norm(d);
if dnorm == 0
    error('linesch_ww_mod: d is zero')
end
t = 1;  % important to try steplength one first
nfeval = 0;
nbisect = 0;
nexpand = 0;
% the following limits are rather arbitrary
% nbisectmax = 30; % 50 is TOO BIG, because of rounding errors
nbisectmax = max(30, round(log2(1e5*dnorm))); % allows more if ||d|| big
nexpandmax = max(10, round(log2(1e5/dnorm))); % allows more if ||d|| small
done = 0;
while ~done
    x = x0 + t*d;
    nfeval = nfeval + 1;
    [f,grad] = feval(fgname, x, pars);
    fevalrec(nfeval) = f;
    if f < fvalquit % nothing more to do, quit
        fail = 0;
        alpha = t; % normally beta is inf
        xalpha = x;
        falpha = f;
        gradalpha = grad;
        return
    end
    gtd = grad'*d;
    % the first condition must be checked first. NOTE THE >=.
    if f >= f0 + c1*t*g0 | isnan(f) % first condition violated, gone too far
        beta = t; 
        gradbeta = grad; % discard f
    % now the second condition.  NOTE THE <=
    elseif gtd <= c2*g0 | isnan(gtd) % second condition violated, not gone far enough
        alpha = t;
        xalpha = x;
        falpha = f;
        gradalpha = grad;
    else   % quit, both conditions are satisfied
        fail = 0;
        alpha = t;
        xalpha = x;
        falpha = f;
        gradalpha = grad;
        beta = t;
        gradbeta = grad;
        return
    end
    % setup next function evaluation
    if beta < inf
        if nbisect < nbisectmax
            nbisect = nbisect + 1;
            t = (alpha + beta)/2; % bisection
        else
            done = 1;
        end
    else
        if nexpand < nexpandmax
            nexpand = nexpand + 1;
            t = 2*alpha;  % still in expansion mode
        else
            done = 1;
        end
    end 
end % loop
% Wolfe conditions not satisfied: there are two cases
if beta == inf % minimizer never bracketed
    fail = -1;
    if prtlevel > 1
        fprintf('Line search failed to bracket point satisfying weak ');
        fprintf('Wolfe conditions, function may be unbounded below\n')
    end
else % point satisfying Wolfe conditions was bracketed
    fail = 1;
    if prtlevel > 1
        fprintf('Line search failed to satisfy weak Wolfe conditions')
        fprintf(' although point satisfying conditions was bracketed\n')
    end
end

function ini = isnaninf(M)
% returns the scalar 1 if ANY entry of M is nan or inf; 0 otherwise
% note: isnan and isinf return matrices if M is a matrix, and
% if treats [0 1] as false, not true.
% ini = sum(sum(isnan(M))) > 0 | sum(sum(isinf(M))) > 0;
ini = any(any(isnan(M))) | any(any(isinf(M)));
%
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = hgprod(H0, g, S, Y)
%  compute the product required by the LM-BFGS method
%  see Nocedal and Wright
%   Send comments/bug reports to Michael Overton, overton@cs.nyu.edu,
%   with a subject header containing the string "hanso" or "gradsamp".
%   Version 2.0, 2010, see GPL license info below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(S,2);  % number of saved vector pairs (s,y) 
q = g;
for i = N:-1:1
   s = S(:,i);
   y = Y(:,i);
   rho(i) = 1/(s'*y);
   alpha(i) = rho(i)*(s'*q);
   q = q - alpha(i)*y;
end
r = H0*q;
for i=1:N
   s = S(:,i);
   y = Y(:,i);
   beta = rho(i)*(y'*r);
   r = r + (alpha(i)-beta)*s;
end

function [x,d,q,info] = qpspecial(G,varargin)
% Call:
% [x,d,q,info] = qpspecial(G,varargin)
%
% Solves the QP
%
%    min     q(x)  = || G*x ||_2^2 = x'*(G'*G)*x
%    s.t.  sum(x)  = 1
%              x  >= 0
%
% The problem corresponds to finding the smallest vector
% (2-norm) in the convex hull of the columns of G
%
% Inputs:
%     G            -- (M x n double) matrix G, see problem above
%     varargin{1}  -- (int) maximum number of iterates allowed
%                     If not present, maxit = 100 is used
%     varargin{2}  -- (n x 1 double) vector x0 with initial (FEASIBLE) iterate.
%                     If not present, (or requirements on x0 not met) a
%                     useable default x0 will be used
%
% Outputs:
%     x       -- Optimal point attaining optimal value
%     d = G*x -- Smallest vector in the convex hull
%     q       -- Optimal value found = d'*d
%     info    -- Run data:
%                info(1) =
%                   0 = everything went well, q is optimal
%                   1 = maxit reached and final x is feasible. so q
%                       might not be optimal, but it is better than q(x0)
%                   2 = something went wrong
%                info(2) = #iterations used
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 2.0 Copyright (C) 2010  Anders Skajaa
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo  = 0;                    % set echo = 1 for printing 
                              % (for debugging). Otherwise 0. 
                              %
[m,n] = size(G);              % size of problem
if ~(m*n>0)                   % in this case
    fprintf(['qpspecial:',... % G is empty, so nothing we can do
        ' G is empty.\n']);   % exit with warning 
    info = [2,0];             % info(1) = 2;
    x = []; d = []; q = inf;  % and empty structs
    return;                   % and optimal value is inf
end                           %
                              %
e     = ones(n,1);            % vector of ones
                              %
if nargin > 1                 % set defauls
    maxit = varargin{1};      % maximal # iterations
    maxit = ceil(maxit);      % in case of a non-integer input
    maxit = max(maxit,5);     % always allow at least 5 iterations
else                          %
    maxit = 100;              % default is 100
end                           % which is always plenty
if nargin > 2                 % 
    x  = varargin{2};         % if x0 is specified
    x  = x(:);                % if given as row instead of col
    nx = size(x,1);           % check that size is right
    if any(x<0) || nx~=n      % use it, unless it is 
        x = e;                % infeasible in the ineq 
    end                       % constraints. 
else                          % use it, otherwise
    x = e;                    % use (1,1,...,1)
end                           % which is an interior point
                              %
idx   = (1:(n+1):(n^2))';     % needed many times
Q     = G'*G;                 % Hessian in QP
z     = x;                    % intialize z
y     = 0;                    % intialize y
eta   = 0.9995;               % step size dampening
delta = 3;                    % for the sigma heuristic
mu0   = (x'*z) / n;           % first my
tolmu = 1e-5;                 % relative stopping tolerance, mu
tolrs = 1e-5;                 % and residual norm
kmu   = tolmu*mu0;            % constant for stopping, mu
nQ    = norm(Q,inf)+2;        % norm of [Q,I,e]
krs   = tolrs*nQ;             % constant for stopping, residuals
ap    = 0; ad = 0;            % init steps just for printing
if echo > 0                   % print first line
    fprintf(['k    mu   ',... %
   '     stpsz     res\n',... %
   '-----------------',...    %
   '-----------------\n']);   % 
end                           %
                              %
for k = 1:maxit               %
                              %
    r1 = -Q*x + e*y + z;      % residual
    r2 = -1 + sum(x);         % residual
    r3 = -x.*z;               % slacks
    rs = norm([r1;r2],inf);   % residual norm
    mu = -sum(r3)/n;          % current mu
                              %
                              % printing (for debugging)
    if echo > 0
        fprintf('%-3.1i %9.2e %9.2e %9.2e \n',...
            k,mu/mu0,max(ap,ad),rs/nQ);
    end
                              % stopping
    if mu < kmu               % mu must be small
        if rs < krs           % primal feas res must be small
            info = [0,k-1];   % in this case, all went well
            break;            % so exit with info = 0
        end                   % so exit loop
    end                       %
                              %
    zdx     = z./x;           % factorization
    QD      = Q;              %  
    QD(idx) = QD(idx) + zdx;  % 
    [C,ef]  = chol(QD);       % C'*C = QD
    if ef > 0                 % safety to catch possible
        info = [2,k];         % problems in the choleschy  
        break;                % in this case, 
    end                       % break with info = 2
    KT      = C'\e;           % K'   = (C')^(-1) * e
    M       = KT'*KT;         % M    = K*K'
                              %
    r4  = r1+r3./x;           % compute approx 
    r5  = KT'*(C'\r4);        % tangent direction
    r6  = r2+r5;              % using factorization 
    dy  = -r6/M;              % from above
    r7  = r4 + e*dy;          %
    dx  = C\(C'\r7);          %
    dz  = (r3 - z.*dx)./x;    %
                              %
    p   = -x ./ dx;           % Determine maximal step 
    ap  = min(min(p(p>0)),1); % possible in the 
    if isempty(ap)            % approx tangent direction
        ap = 1;               % here primal step size
    end                       % 
    p   = -z ./ dz;           % here dual step size
    ad  = min(min(p(p>0)),1); % 
    if isempty(ad)            % using different step sizes
        ad = 1;               % in primal and dual improves
    end                       % performance a bit
                              %  
    muaff = ((x + ap*dx)'*... % heuristic for 
             (z + ad*dz))/n;  % the centering parameter
    sig   = (muaff/mu)^delta; %
                              %
    r3  = r3 + sig*mu;        % compute the new corrected
    r3  = r3 - dx.*dz;        % search direction that now 
    r4  = r1+r3./x;           % includes the appropriate
    r5  = KT'*(C'\r4);        % amount of centering and 
    r6  = r2+r5;              % mehrotras second order 
    dy  = -r6/M;              % correction term (see r3).
    r7  = r4 + e*dy;          % we of course reuse the 
    dx  = C\(C'\r7);          % factorization from above
    dz  = (r3 - z.*dx)./x;    %
                              %
    p   = -x ./ dx;           % Determine maximal step 
    ap  = min(min(p(p>0)),1); % possible in the 
    if isempty(ap)            % new direction
        ap = 1;               % here primal step size
    end                       % 
    p   = -z ./ dz;           % here dual step size
    ad  = min(min(p(p>0)),1); % 
    if isempty(ad)            % 
        ad = 1;               % 
    end                       % 
                              % update variables
    x   = x + eta * ap * dx;  % primal
    y   = y + eta * ad * dy;  % dual multipliers
    z   = z + eta * ad * dz;  % dual slacks 
                              %
end                           % end main loop
                              %
if k == maxit                 % if reached maxit
    info = [1,k];             % set info(1) = 1
end                           %
x = max(x,0);                 % project x onto R+
x = x/sum(x);                 % so that sum(x) = 1 exactly                        
d = G*x;                      % and set other output 
q = d'*d;                     % variables using best 
                              % found x
                              %
if echo > 0                   % printing
    str = 'optimal';          % status string to print
    if info(1)==1             % in the last line
        str = 'maxit reached';%
    elseif info(1)==2         %
        str = 'failed';       %
    end                       %
    fprintf(['---------',...  % print last line
    '-------------------',... %
    '------\n',...            %     
    ' result: %s \n',...      %
    '-------------------',... %
    '---------------\n'],str);%
end                           %
    




