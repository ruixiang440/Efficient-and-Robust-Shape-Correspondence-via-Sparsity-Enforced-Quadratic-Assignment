function [x, g, out]= fminPGBBgpu(x, fun, proj, opts, varargin)
%-------------------------------------------------------------------------
% Line search algorithm for optimization on manifold:
%
%   min f(X), where X \in R^{n,p} and X \in U
%
%
% Input:
%           X --- 
%         fun --- objective function and its gradient:
%                 [F, G] = fun(X,  data1, data2)
%                 F, G are the objective function value and gradient, repectively
%                 data1, data2 are addtional data, and can be more
%                 Calling syntax:
%                   [X, out]= OptManiMulitBallGBB(X0, @fun, opts, data1, data2);
%
%        proj --- projection operator
%
%        opts --- option structure with fields:
%                 record = 0, no print out
%                 mxiter       max number of iterations
%                 xtol        stop control for ||X_k - X_{k-1}||
%                 gtol        stop control for the projected gradient
%                 ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
%                             usually, max{xtol, gtol} > ftol
%   
% Output:
%           x --- solution
%           g --- gradient of x
%         Out --- output information
%
%
% Author: Zaiwen Wen
%   Version 1.0 .... 2010/10
%-------------------------------------------------------------------------


if nargin < 2; error('call [x, g, out]= fminPGBB(x, fun, proj, opts)'); end
if nargin < 3; proj = @(x)x; end
if nargin < 4; opts = [];    end

if isempty(proj); proj = @(x)x; end
   
% termination rule
if ~isfield(opts, 'xtol');      opts.xtol = 1e-10; end
if ~isfield(opts, 'gtol');      opts.gtol = 1e-6; end
if ~isfield(opts, 'ftol');      opts.ftol = 1e-14; end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'tau');       opts.tau  = 1e-3; end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-4; end
if ~isfield(opts, 'eta');       opts.eta  = 0.1; end
if ~isfield(opts, 'gamma');     opts.gamma  = 0.85; end
if ~isfield(opts, 'STPEPS');    opts.STPEPS  = 1e-10; end
if ~isfield(opts, 'nt');        opts.nt  = 5; end
if ~isfield(opts, 'maxit');     opts.maxit  = 2000; end
if ~isfield(opts, 'eps');       opts.eps = 1e-6; end
if ~isfield(opts, 'record');    opts.record = 0; end

%-------------------------------------------------------------------------------
% copy parameters
xtol = opts.xtol;
ftol = opts.ftol;
gtol = opts.gtol;
maxit = opts.maxit;
rhols = opts.rhols;
eta   = opts.eta;
gamma = opts.gamma;
record = opts.record;
nt = opts.nt;
% crit = gpuArray(ones(nt, 3)); %%%%%%%%%%
crit = ones(nt, 3);


[n,p] = size(x);

%% Initial function value and gradient
% prepare for iterations
[f,g] = feval(fun, x, varargin{:});    out.nfe = 1;  out.fval0 = f;
nrmG  = norm(full(g), 'fro');

Q = 1; Cval = f; tau = opts.tau;
%% Print iteration header if debug == 1
if (record >= 1)
    fprintf('----------- fminBB ----------- \n');
    fprintf('%4s \t %10s \t %10s \t  %10s \t %5s \t %9s \t %7s \n', ...
        'Iter', 'tau', 'f(X)', 'nrmG', 'Exit', 'funcCount', 'ls-Iter');
end

if record == 10; out.fvec = f; end
out.msg = 'exceed max iteration';

xbest = feval(proj, x); 



fbest = f;  gbest = g;   nrmGbest = nrmG;
nrmy = 1; sy = 1;

%% main iteration
for iter = 1 : maxit
    xp = x;     fp = f;     gp = g;   
    nls = 1; deriv = rhols*nrmG^2;
    while 1
        % calculate g, f,
        x = feval(proj, xp - tau*gp);
%         x = sparsify(full(x));
        out.tau(iter)=gather(tau);
        %f = feval(fun, x, varargin{:});   out.nfe = out.nfe + 1;
        [f,g] = feval(fun, x, varargin{:});   out.nfe = out.nfe + 1;
        out.e(iter)=f;
        if f <= Cval - tau*deriv || nls >= 5
            break
        end
        tau = eta*tau;
        nls = nls+1;
    end  
    
     if record == 10; out.fvec = [out.fvec; f]; end
    
    if f < fbest
        xbest = x; fbest = f;  gbest = g;
    end
    
    %nrmG  = norm(g, 'fro'); 
    nrmG  = norm(x - feval(proj, x - g), 'fro'); 
    s = x - xp; XDiff = norm(s,'fro')/sqrt(n);
    FDiff = abs(fp-f)/(abs(fp)+1);

    if (record >= 1)
        fprintf('%4d  %3.2e  %7.6e  %3.2e  %3.2e  %3.2e  %2d\n', ...
            iter, tau, f, nrmG, XDiff, FDiff, nls);
    end
    
    crit(iter,:) = [gather(nrmG), gather(XDiff), FDiff];
    mcrit = mean(crit(iter-min(nt,iter)+1:iter, :),1);
    
    %if nrmG < gtol
    if ( XDiff < xtol && FDiff < ftol ) || nrmG < gtol ...
           % || all(mcrit(2:3) < 10*[xtol, ftol])  
        out.msg = 'converge';
        break;
    end
    
    y = g - gp;   nrmy = norm(y,'fro');    
    sy = abs(iprod(s,y));    tau = opts.tau;
%     if sy > 1e-16 && nrmy >= 1e-12 ; 
        if mod(iter,2)==0; tau = norm(s,'fro')^2/sy;
        else tau = sy/(nrmy^2); end
%     else
%         tau = 0.001*norm(x,'fro')/nrmG;
%     end
    % safeguarding on tau
    tau = max(min(tau, 1e20), 1e-20);
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + f)/Q;
end

x = xbest;  g = gbest;
out.nrmG = nrmGbest;
out.fval = fbest;
out.iter = iter;
end

function a = iprod(x,y)
    a = full(sum(sum(x.*y)));
end

