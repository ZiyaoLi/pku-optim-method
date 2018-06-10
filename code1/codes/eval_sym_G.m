function [ G ] = eval_sym_G( x,symG )
% EVAL_SYM_G.M function evaluates matrix function value G(x) of a given 
% function symbolic expression on the point x. It can be a lot slower than 
% directly calculating the function.
%
% Input:
% x:        x in f, n-vector. Note that x must be column vectors.
% symG:     Given function expression of G, sym.
%           The expression shall calculate G(x) as G = G(x1, x2, ..., xn),
%           where x_i's are symbols, too. Note that given x must have the 
%           same dimension as the variables in symf do.% 
% Output:
% G:        The calculated function value G(x), matrix.
% 
% Call:
% [G] = eval_diff_G(x, symG)

% Date:     Apr 11th, 2018
% Creator:  BroC


n = length(x);
X = sym('x',[1,n]);
G = eval(subs(symG,X,x'));

end

