function [ x, info, perf ] = ...
    global_bb( fg_handle, Args, x0, Std, Rule, varargin )
% CONJUGATE_GRADIENT.M finds the answer to problem P: x* = argmin f(x) 
% using conjugate gradient method. 
% 
% Input 
% fg_handle:	f & g in P, function_handle.
%               The function should be declared as a M-function
%                   [F g] = fg_handle(x; p1, p2, ..., pt),
%               where f is a scalar, along with g an n-vector.
% Args:         arguments for the global BB method, 7-vector
%       Args(1):    alpha_0, default=1.
%       Args(2):    M, positive integer, default=5.
%       Args(3):    gamma in (0, 1), default=0.8.
%       Args(4):    delta, default=1.
%       Args(5):    sigma, 0 < sigma < 1, default=0.5.
%       Args(6):    eps, default=1e-8.
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

narginchk(3, inf);
nargoutchk(0, 3);

% setting defaults.
if nargin == 2 || isempty(Std)
    Std = [1e-8, 1e-6, 100, inf];
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

alpha = Args(1);
M = Args(2);
gamma = Args(3);
delta = Args(4);
sigma = Args(5);
eps = Args(6);
Eps = 1 / eps;

delta_f = 0;
tp = g' * g;
g_norm = norm(g, pnorm);
% c = g_norm * eps2;

f_hist = ones(M, 1) .* (-Inf);
g_hist = ones(M, 1) .* Inf;

while ( delta_f > eps1 || g_norm > eps2 ) && iter < maxiter
    iter = iter + 1;
    if verbose && rem(iter,50) == 0
        fprintf('Iteration %d: delta=%3e, g_norm=%3e, f=%3e\n',...
            iter, delta_f, g_norm, f);
    end
    if record_fg
        perf.f_rec(iter+1) = f;
        perf.g_norm_rec(iter+1) = norm(g, Std(4));
    end
    
    d = -g;
    
    % circulate storage to update M-history of f and |g|_2 value.
    f_hist(rem(iter - 1, M) + 1) = f;
    g_hist(rem(iter - 1, M) + 1) = gamma * tp;
    
    if alpha <= eps || alpha >= Eps
        alpha = delta;
    end
    
    while feval(fg_handle, x + alpha .* d, false, varargin{:}) ...
            > max(f_hist - alpha .* g_hist)  % S4 -> S6
        alpha = sigma * alpha;
        count = count + 1;
    end
      
    % break wrt. wrong step length.
%     if ls_info(1) == -1
%         fprintf(['Error: No step can be found on this direction. '...
%             'CG method cannot continue at iteration %d.\n'], iter)
%         info(1) = -1;
%         break
%     end
    
    x = x + alpha .* d;  % S5
    % norm(x, Inf)
    
    fp = f; % to remember the previous values
    gp = g;
    
    [f, g] = feval(fg_handle, x, true, varargin{:});
    alpha =  (alpha * tp) / (tp - gp'*g);
    
    tp = g' * g; % update the storaged ||g||_2
    delta_f = abs(f - fp);
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
info(3) = count + 2 * iter; % reducing count++ in the loop by directly
                        % adding iter to the counter, since 2 eval per 
                        % iter - remind the one escaping from the inner 
                        % loop which did not call the counter.
perf.x = x;
perf.f = f;
perf.g = g;

end