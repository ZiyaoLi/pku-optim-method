function [ x, info, perf ] = ...
    lm_newton( fg_handle, G_handle, x0, Std, Rule, varargin )
% LM_NEWTON.M finds the answer to problem P: x* = argmin f(x) using
% LM fixed newton method. 
% 
% Input 
% fg_handle:	f & g in P, function_handle.
%               The function should be declared as a M-function
%                   [F g] = fg_handle(x; p1, p2, ..., pt),
%               where f is a scalar, along with g an n-vector.
% G_handle:     G in P, function_handle.
%               The function should be declared as a M-function
%                   [G] = G_handle(x, q1, q2, ...),
%               where G is a n x n matrix, the Hesse matrix of function f.
% x0:           x0 in P, the initial value of x, n-vector
% Std:         4-vector
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
%               0 - Successful call.
%               1 - Reach the maximum of iterations.
%              -1 - Method cannot continue.
%       info(2):    the numbers of iterations.
%       info(3):    the numbers of evaluating function fg_handle.
% perf:     other useful data, struct.
%       perf.x:     new point calculated with the method
%       perf.f:     function value at the point perf.x.
%       perf.g:     gradient at the point perf.x.
% 
% Call
% [x,info,perf] = lm_newton(fg_handle, G_handle, x0)
% [x,info,perf] = lm_newton(fg_handle, G_handle, x0, Std)
% [x,info,perf] = lm_newton(fg_handle, G_handle, x0, Std, Rule)
% [x,info,perf] = lm_newton(fg_handle, G_handle, x0, Std, Rule, ...
%                     n_param, p1, p2, ..., p_{n_param}, q1, q2, ...)

% Date:     Apr 10th, 2018
% Creator:  broC


narginchk(3, inf);
nargoutchk(0, 3);

% setting defaults.
if nargin == 3 || isempty(Std)
    Std = [1e-8, 1e-8, 100, inf];
end
if nargin <= 4 || isempty(Rule)
    Rule.opt = 1;
    Rule.opt = bodfltchk(Rule.opt, [1 10 10 0.95 0.05]);
    Rule.crtr = bodfltfunc(Rule.crtr, @bostwlf);
    Rule.mthd = bodfltfunc(Rule.mthd, @bointrplt33);
end

% consider the matrix H = G + aI, where H is a positive definite matrix and
% G is not. G has its max eigenvalue > 0 and minimal one < 0. Find scalar a
% s.t. cond(H) = max_cond so that the problem is not ill-defined. The
% method is more robust than simply adding 1 on the G diagonal.
max_cond = 20;

n = length(x0);
info = [0 0 0];
n_param_fun = varargin{1}+1;
n_param_total = length(varargin);
if n_param_fun > n_param_total
    error('Wrong parameter number provided.')
end
x = x0;
i = 0;
[f,g] = feval(fg_handle, x, varargin{2:n_param_fun});
n = length(x);
delta = Std(1) + 1; % to avoid ending before the first iteration.
g_norm = norm(g, Std(4));

while ( delta > Std(1) || g_norm > Std(2) ) && i < Std(3)
    i = i + 1;
    if rem(i,10) == 0
        fprintf('Iteration %d\n',i);
    end
    perf.f_rec(i) = f;
    perf.g_norm_rec(i) = norm(g, Std(4));
    delta = f; % only to remember the previous f value
    G = feval(G_handle, x, varargin{n_param_fun+1:n_param_total});
    e = eig(G);
    if min(e) <= 0
        fprintf('G not SPD at iteration %d\n',i)
        c = (max(e) - max_cond * min(e)) / (max_cond - 1);
        G = G + c .* eye(n);
    end
    d = - G \ g;
    
    [alpha, rst] = bolinesearch(fg_handle, x, d, Rule, varargin{2:n_param_fun});
    info(3) = info(3) + rst(3);
    if rst(1) == -1
        fprintf(['Error: No step can be found on this direction. '...
            'Newton method cannot continue at iteration %d.\n'], i)
        info(1) = -1;
        info(2) = i;
        return
    end
    x = x + alpha .* d;
    [f,g] = feval(fg_handle, x, varargin{2:n_param_fun});
    delta = abs(f - delta);
    g_norm = norm(g, Std(4));
end

% output relevant information
perf.f_rec(i+1) = f;
perf.g_norm_rec(i+1) = norm(g, Std(4));
info(1) = 0;
info(2) = i;
info(3) = info(3) + (i + 1) * (n + 1);
perf.x = x;
perf.f = f;
perf.g = g;
if delta > Std(1) || g_norm > Std(2)
    % indicating that the iteration ended because the max_iter reached.
    info(1) = 1;
end

end