function [ x, info, perf ] = ...
    conjugate_gradient( fg_handle, Method, x0, Std, Rule, varargin )
% CONJUGATE_GRADIENT.M finds the answer to problem P: x* = argmin f(x) 
% using conjugate gradient method. 
% 
% Input 
% fg_handle:	f & g in P, function_handle.
%               The function should be declared as a M-function
%                   [F g] = fg_handle(x; p1, p2, ..., pt),
%               where f is a scalar, along with g an n-vector.
% Method:       indicator of the CG formula used, scalar.
%               1 - FR    (Fletcher-Reeves)
%               2 - PRP   (Polak-Ribiere-Polyak)
%               3 - PRP+  (Polak-Ribiere-Polyak Plus)
%               4 - CD    (Conjugate Descend)
%               5 - DY    (Dai-Yuan)
% x0:           x0 in P, the initial value of x, n-vector
% Std:          4-vector
%       Std(1):     the epsilon of the stop criterion |f_k-f_{k-1}| < eps.
%                   default 1e-8.
%       Std(2):     the epsilon of the stop criterion: ||g_k|| < eps.
%                   default 1e-8.
%       Std(3):     the maximum iteration allowed in this program.
%                   default 100.
%       Std(4):     the p-norm used in the stop criterion ||g_k||_p < eps.
%                   typically p = 1, 2, inf. default inf.
% Rule:     option & method & criterion to solve P, struct.
%           To see the detailed description, please check BOLINESEARCH.M 
% n_param:  the number of parameters used in fg_handle, i.e. the t in the 
%           function form f(x; p1, p2, ..., pt), scalar.
%
% Output
% x:        x* in P, the optimal answer, n-vector.
% info:     status to execute this function, 3-vector.
%       info(1):    exit code.
%                   0 - Successful call.
%                   1 - Reach the maximum of iterations.
%                  -1 - Method cannot continue.
%       info(2):    the numbers of iterations.
%       info(3):    the numbers of evaluating function fg_handle.
% perf:     other useful data, struct.
%       perf.x:     new point calculated with the method
%       perf.f:     function value at the point perf.x.
%       perf.g:     gradient at the point perf.x.
% 
% Call
% [x,info,perf] = quasi_newton(fg_handle, Method, x0)
% [x,info,perf] = quasi_newton(fg_handle, Method, x0, Std)
% [x,info,perf] = quasi_newton(fg_handle, Method, x0, Std, Rule)
% [x,info,perf] = quasi_newton(fg_handle, Method, x0, Std, Rule, ...
%                     n_param, p1, p2, ..., p_{n_param}, q1, q2, ...)

% Date:     Jun 6th, 2018
% Creator:  broC

% control whether detailed iterations are shown; whether the variation of f
% and g_norm are monitored.
verbose = true;
record_fg = false;
n_restart = 5;

narginchk(3, inf);
nargoutchk(0, 3);

% setting defaults.
if nargin == 2 || isempty(Std)
    Std = [1e-8, 1e-8, 100, inf];
end
if nargin <= 4 || isempty(Rule)
    Rule.opt = 1;
    Rule.opt = bodfltchk(Rule.opt, [1 10 10 0.95 0.05]);
    Rule.crtr = bodfltfunc(Rule.crtr, @bostwlf);
    Rule.mthd = bodfltfunc(Rule.mthd, @bointrplt33);
end

info = [0 0 0];

x = x0;
iter = 0;
[f,g] = feval(fg_handle, x, true, varargin{:});
count = 1;

eps1 = Std(1);
eps2 = Std(2);
maxiter = Std(3);
pnorm = Std(4);

delta = 0;
g_norm = norm(g, pnorm);
% c = g_norm * eps2;

while ( delta > eps1 || g_norm > eps2 ) && iter < maxiter
    iter = iter + 1;
    if verbose && rem(iter,50) == 0
        fprintf('Iteration %d: delta=%3e, g_norm=%3e, f=%3e\n',...
            iter, delta, g_norm, f);
    end
    
    if record_fg
        perf.f_rec(iter+1) = f;
        perf.g_norm_rec(iter+1) = norm(g, Std(4));
    end
    
    if rem(iter, n_restart) == 1 % iter == 1%
        d = -g;
        tp = g' * g;
    else
        switch(Method)
            % record last-time calculation in tp
            case 1
                t = g' * g;  beta = t / tp;  tp = t;
            case 2
                t = g' * g;  beta = (t - g' * gp) / tp;  tp = t;
            case 3
                t = g' * g;  beta = max((t - g' * gp) / tp, 0);  tp = t;
            case 4
                beta = (g' * g) / (d' * gp);
            case 5
                beta = (g' * g) / (d' * (g - gp));
        end
        d = -g + beta .* d;
    end
    
    [alpha, ls_info] = bolinesearch(fg_handle, x, d, Rule, varargin{:});
    count = count + ls_info(3);
    
    % break wrt. wrong step length.
    if ls_info(1) == -1
        fprintf(['Error: No step can be found on this direction. '...
            'CG method cannot continue at iteration %d.\n'], iter)
        info(1) = -1;
        break
    end
    
    x = x + alpha .* d;
    % norm(x, Inf)
    fp = f; % only to remember the previous f value
    gp = g;
    [f, g] = feval(fg_handle, x, true, varargin{:});
    delta = abs(f - fp);
    g_norm = norm(g, pnorm);
end

% output relevant information
if record_fg
    perf.f_rec(iter+1) = f;
    perf.g_norm_rec(iter+1) = norm(g, Std(4));
end
if iter == maxiter
    info(1) = 1;
end

info(2) = iter;
info(3) = count + iter; % reducing count++ in the loop by directly adding
                        % iter to the counter, since one eval per iter.
perf.x = x;
perf.f = f;
perf.g = g;

end