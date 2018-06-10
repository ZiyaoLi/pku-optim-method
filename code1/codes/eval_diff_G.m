function [ G ] = eval_diff_G( x, fg_handle, eps, varargin )
% EVAL_DIFF_FG.M function calculates function Hesse matrix G(x) of a given 
% function & gradient handle on the point x using differential method.
%
% Input:
% x:        x in f, n-vector.
% fg_handle:    Given function & gradient handle of f, function_handle.
%           The function & gradient handle should be arranged as
%               [f g] = fg_handle(x; p1, p2, ...).
% eps:      The differential step of the differential method, scalar, i.e.
%                   g(x_i) <- (f(x + eps .* e_i) - f(x)) / eps,
%           default 1e-8.
% varargin: Other parameters, p1, p2, ..., to pass to function f. i.e.
%                   f = f(x; p1, p2, ...).
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
nargoutchk(0, 3);
% setting defaults.
if nargin == 2 || isempty(eps)
    eps = 1e-8;
end

n = length(x);
[~, g] = feval(fg_handle, x, varargin{:});
G = zeros(n);
for i=1:n
    xi = x;
    xi(i) = x(i) + eps;
    [~,gi] = feval(fg_handle, xi, varargin{:});
    G(i,1:n) = (gi - g) ./ eps;
end
end

