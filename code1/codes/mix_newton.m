function [ x, info, perf ] = ...
    mix_newton( ObjFun, Hesse, Point, Std, Rule, Params, varargin )
%NEWTON Find the answer to problem P: x* = argmin f(x) using basic newton
% method. That is alpha = 1 for iteration x <- x + alpha * d
%
% Call
% [ xStar, info, perf ] = newton( ObjFun, Hesse, Point )
% [ xStar, info, perf ] = newton( ObjFun, Hesse, Point, Std )
% [ xStar, info, perf ] = newton( ObjFun, Hesse, Point, Std, Rule )
% [ xStar, info, perf ] = newton( ObjFun, Hesse, Point, Std, Rule,...
%                                 n_param_f, p1, p2, ..., q1, q2, ... )
%
% Input 
% ObjFun:   f & g in P, function_handle
%               The function should be declared as a M-function
%                   [F g] = f(x, p1,p2,...)
%               F is a scalar, along with g an n-vector.
% Hesse:    G in P, function_handle
%               The function should be declared as a M-function
%                   [G] = f(x, q1,q2,...)
%               F is a scalar, along with g an n-vector.
% Point:    x in P, n-vector
% Step:     4-vector
%       Step(1):    the epsilon of the stop criterion |f_k-f_{k-1}| < eps.
%                   default 1e-8.
%       Step(2):    the epsilon of the stop criterion: ||g_k|| < eps.
%                   default 1e-8.
%       Step(3):    the maximum iteration allowed in this program.
%                   default 100.
%       Step(4):    the p-norm used in the stop criterion ||g_k||_p < eps.
%                   typically p = 1, 2, inf. default inf.
% Rule:     option & method & criterion to solve P, struct
%       To see the detailed description, please check bolinesearch.m 
% n_param:  the number of parameters used in ObjFun, scalar
%
% Output
% x:    x in P the optimal answer, n-vector
% info:     status to execute this function, integral array
%       info(1):    exit code
%               0 - Successful call
%               1 - Reach the maximum of iterations
%              -1 - a is not real
%       info(2):    the recorded F value during the iteration.
%       info(3):    the recorded g norm during the iteration.
% perf:     other useful data, struct
%       perf.x:     new point after iteration
%       perf.F:     function value at new point
%       perf.g:     gradient at new point

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

info = [0 0 0];
n_param_fun = varargin{1}+1;
n_param_total = length(varargin);
if n_param_fun > n_param_total
    error('Wrong parameter number provided.')
end
x = Point;
i = 0;
[f,g] = feval(ObjFun, x, varargin{2:n_param_fun});
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
    G = feval(Hesse, x, varargin{n_param_fun+1:n_param_total});
    if cond(G) > 1e8
        d = - g;
        perf.case(i) = -1;
    else
        d = - G \ g;
        perf.case(i) = 0;
        if g' * d > Params(1) .* (g' * g) * (d' * d)
            % d is so close to g that the function may not decrease.
            d = -d;
            perf.case(i) = 1;
        elseif abs(g' * d) > Params(2) .* (g' * g) * (d' * d)
            % orthogonal direction with gradient
            d = - g;
            perf.case(i) = 2;
        end
    end
    
    [alpha, rst] = bolinesearch(ObjFun, x, d, Rule, varargin{2:n_param_fun});
    perf.a_rec(i) = alpha;
    info(3) = info(3) + rst(3);
    if rst(1) == -1
        fprintf(['Error: No step can be found on this direction. '...
            'Newton method cannot continue at iteration %d.\n'], i)
        info(1) = -1;
        info(2) = i;
        return
    end
    x = x + alpha .* d;
    [f,g] = feval(ObjFun, x, varargin{2:n_param_fun});
    delta = abs(f - delta);
    g_norm = norm(g, Std(4));
end

% output relevant information
info(1) = 0;
info(2) = i;
info(3) = info(3) + i + 1;
perf.x = x;
perf.F = f;
perf.g = g;
if delta > Std(1) || g_norm > Std(2)
    % indicating that the iteration ended because the max_iter reached.
    info(1) = 1;
end

end