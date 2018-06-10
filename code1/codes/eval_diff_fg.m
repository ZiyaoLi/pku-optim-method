function [ f,g ] = eval_diff_fg( x, f_handle, eps, varargin )
% EVAL_DIFF_FG.M function calculates function value f(x) and gradient g(x)
% of a given function handle on the point x using differential method.
%
% Input:
% x:        x in f, n-vector.
% f_handle: Given function handle of f, function_handle.
% eps:      The differential step of the differential method, scalar, i.e.
%                   g(x_i) <- (f(x + eps .* e_i) - f(x)) / eps,
%           default 1e-8.
% varargin: Other parameters to pass to function f. i.e.
%                   f = f(x; p1, p2, ...)
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = eval_diff_fg(x, f_handle)
% [f, g] = eval_diff_fg(x, f_handle, eps)
% [f, g] = eval_diff_fg(x, f_handle, eps, p1, p2, ...)

% Date:     Apr 11th, 2018
% Creator:  BroC


narginchk(2, inf);
nargoutchk(0, 2);
% setting defaults.
if nargin == 2 || isempty(eps)
    eps = 1e-8;
end

f = feval(f_handle, x, varargin{:});
for i=1:length(x)
    xi = x;
    xi(i) = x(i) + eps;
    g(i) = (feval(f_handle, xi, varargin{:}) - f) / eps;
end
g = g';
end
