function [ f,g ] = eval_sym_fg( x,symf,symg )
% EVAL_SYM_FG.M function evaluates function value f(x) and gradient g(x)
% of a given function symbolic expression on the point x. It can be a lot
% slower than directly calculating the function and gradient value.
%
% Input:
% x:        x in f, n-vector. Note that x must be column vectors.
% symf:     Given function expression of f, sym.
%           The expression shall calculate f(x) as f = f(x1, x2, ..., xn),
%           where x_i's are symbols, too. Note that given x must have the 
%           same dimension as the variables in symf do.
% symg:     Given gradient expression of g, sym.
%           The expression shall calculate g(x) as g = g(x1, x2, ..., xn).
%           where x_i's are symbols, too. Note that given x must have the 
%           same dimension as the variables in symg do.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = eval_diff_fg(x, symf, symg)

% Date:     Apr 11th, 2018
% Creator:  BroC


n = length(x);
X = sym('x', [1,n]);
f = eval(subs(symf,X,x'));
g = eval(subs(symg,X,x'));
end
