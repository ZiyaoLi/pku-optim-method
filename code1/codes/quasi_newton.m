function [ x, info, perf ] = ...
    quasi_newton( fg_handle, Method, x0, Std, Rule,  varargin )
% QUASI_NEWTON.M finds the answer to problem P: x* = argmin f(x) using
% quasi-newton method. 
% 
% Input 
% fg_handle:	f & g in P, function_handle.
%               The function should be declared as a M-function
%                   [F g] = fg_handle(x; p1, p2, ..., pt),
%               where f is a scalar, along with g an n-vector.
% Method:       indicator of the quasi-newton formula used, 2-vector
%       Method(1):  the indicator of the method used.
%                   0 - SR1
%                   1 - DFG
%                   2 - BFGS
%       Method(2):  the indicator of matrix forms in the method.
%                   0 - use matrix in H form
%                   1 - use matrix in B form
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

% Date:     Apr 10th, 2018
% Creator:  broC


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
M = eye(length(x));
i = 0;
[f,g] = feval(fg_handle, x, varargin{:});
delta = Std(1) + 1; % to avoid ending before the first iteration.
g_norm = norm(g, Std(4));

while ( delta > Std(1) || g_norm > Std(2) ) && i < Std(3)
    i = i + 1;
    if rem(i,10) == 0
        fprintf('Iteration %d\n',i);
    end
    perf.f_rec(i) = f;
    perf.g_norm_rec(i) = norm(g, Std(4));
    perf.x = x;
    perf.F = f;
    perf.g = g;
    
    delta = f; % only to remember the previous f value
    
    if Method(2)
        % B method adopted
        d = - M \ g;
    else
        % H method adopted
        d = - M * g;
    end
    
    [alpha, rst] = bolinesearch(fg_handle, x, d, Rule, varargin{:});
    info(3) = info(3) + rst(3);
    if rst(1) == -1
        % wrong step length.
        fprintf(['Error: No step can be found on this direction. '...
            'Newton method cannot continue at iteration %d.\n'], i)
        info(1) = -1;
        info(2) = i;
        return
    end
    s = alpha .* d;
    y = g;
    x = x + s;
    [f,g] = feval(fg_handle, x, varargin{:});
    y = g - y;
    switch Method(1)
        case 0
            % SR1
            if Method(2)
                tmp = y - M * s;
                M = M + tmp * tmp' ./ (tmp' * s);
            else
                tmp = s - M * y;
                M = M + tmp * tmp' ./ (tmp' * y);
            end
        case 1
            % DFG
            if Method(2)
                tmp_Bs = M * s;
                tmp_ys = s' * y;
                M = M + ((1 + (tmp_Bs' * s) / tmp_ys) / tmp_ys) .* y * y' -...
                    (y * tmp_Bs' + tmp_Bs * y') ./ tmp_ys;
            else
                tmp = M * y;
                M = M + s * s' ./ (s' * y) - tmp * tmp' ./ (tmp' * y);
            end
        case 2
            % BFGS
            if Method(2)
                tmp = M * s;
                M = M + y * y' ./ (y' * s) - tmp * tmp' ./ (tmp' * s);
            else
                tmp_Hy = M * y;
                tmp_ys = y' * s;
                M = M + ((1 + (tmp_Hy' * y) / tmp_ys) / tmp_ys) .* s * s' -...
                    (s * tmp_Hy' + tmp_Hy * s') ./ tmp_ys;
            end
    end
    delta = abs(f - delta);
    g_norm = norm(g, Std(4));
end

% output relevant information
perf.f_rec(i+1) = f;
perf.g_norm_rec(i+1) = norm(g, Std(4));
info(1) = 0;
info(2) = i;
info(3) = info(3) + i + 1;
perf.x = x;
perf.f = f;
perf.g = g;

if delta > Std(1) || g_norm > Std(2)
    % indicating that the iteration ended because the max_iter reached.
    info(1) = 1;
end

end