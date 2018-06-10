function [ x, info, perf ] = ...
    gauss_newton( frj_handle, x0, Std, Rule,  varargin )
% GAUSS_NEWTON.M finds the answer to the least square problem 
%     P: x* = argmin f(x) = \sum_{i=1}^{m} r_i(x)
% using gauss-newton method. 
% 
% Input 
% frj_handle:	f, r & g in P, function_handle.
%               The function should be declared as a M-function
%                   [f r J] = fg_handle(x; p1, p2, ..., pt),
%               where f is a scalar, along with r an m-vector, J an (m*n) 
%               Jacobian matrix.
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
%       perf.f_rec: function values at every iteraton.
%       perf.g_norm_rec:    gradient norms at every iteraton.
%       perf.g_ratio:   final gradient norm / initial gradient norm.
% 
% Call
% [x,info,perf] = gauss_newton(frj_handle, x0)
% [x,info,perf] = gauss_newton(frj_handle, x0, Std)
% [x,info,perf] = gauss_newton(frj_handle, x0, Std, Rule)
% [x,info,perf] = gauss_newton(frj_handle, x0, Std, Rule, p1, p2, ...)

% Date:     Jun 8th, 2018
% Creator:  broC


% control whether detailed iterations are shown; whether the variation of f
% and g_norm are monitored.
verbose = true;
record_fg = true;
n_verbose = 100;

narginchk(2, inf);
nargoutchk(0, 3);

% setting defaults.
if nargin == 2 || isempty(Std)
    Std = [1e-8, 1e-8, 100, inf];
end
if nargin <= 3 || isempty(Rule)
    Rule.opt = [0, 100, 1, 0.1, 2, 1e-8];
end

info = [0 0 0];

x = x0;
iter = 0;
[f, r, J] = feval(frj_handle, x, true, varargin{:});
count = 1;

eps1 = Std(1);
eps2 = Std(2);
maxiter = Std(3);
pnorm = Std(4);

delta = Std(1) + 1; % to avoid ending before the first iteration.
g = J' * r;
g_norm = norm(g, pnorm);
% c = g_norm * eps2;

while ( delta > eps1 || g_norm > eps2) && iter < maxiter%
    iter = iter + 1;
    if verbose && rem(iter, n_verbose) == 0
        fprintf('Iteration %d: delta=%3e, g_norm=%3e, f=%3e\n',...
            iter, delta, norm(g, pnorm), f);
    end
    
    if record_fg
        perf.f_rec(iter) = f;
        perf.g_norm_rec(iter) = norm(g, Std(4));
    end
    
    d = J \ (-r);
    if any(isnan(d))
        fprintf(['Warning: NaN in direction. '...
            'Taking negative gradient as direction at iteration %d.\n'], iter)
        d = -g;
    end
    
    [alpha, rst] = bolinesearch(frj_handle, x, d, Rule, varargin{:});
    count = count + rst(3);
    if rst(1) == -1
        alpha = 1;
        fprintf(['Warning: No step can be found on this direction. '...
            'Take alpha = 1 at iteration %d.\n'], iter)
    end
    
    x = x + alpha .* d;

    fp = f; % only to remember the previous f value
    [f, r, J] = feval(frj_handle, x, true, varargin{:});
    g = J' * r;
    if any(isnan(g)) || isnan(f) || any(isinf(f)) || any(isinf(r))
        fprintf(['Error: NaN in f or g. '...
            'Gauss-Newton stops at iteration %d.\n'], iter)
        info(1) = -1;
        break
    end
    delta = abs(f - fp);
    g_norm = norm(g, pnorm);
end

% output relevant information
if record_fg
    perf.f_rec(iter+1) = f;
    t = norm(g, pnorm);
    perf.g_norm_rec(iter+1) = t;
    perf.g_ratio = t / perf.g_norm_rec(1);
end
if iter == maxiter
    info(1) = 1;
end

info(2) = iter;
info(3) = count + iter; % reducing count++ in the loop by directly adding
                        % iter to the counter, since one eval per iter.
perf.x = x;
perf.f = f;
perf.r = r;

end