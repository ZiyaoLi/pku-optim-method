% MINIMAL_SURFACE.M solves mininal surface with boundary constraints
% problem with differential method. Detailed description is in "p153",
% Problem 4.3. In this script, the convex set D is limited to a
% square (rectangular) grid, i.e. (0,1)*(0,1) \in R^2.
% 
% Arguments
% M,N:      The number of intervals on one of the differential boundary,
%           i.e. the granularity of the grid, scalar.
% type:     The pre-defined type of boundaries. Please see 'u_bound.m' for
%           detailed description, scalar in {0, 1, 2, 3, 4, 5, 7}.
% V:        The actual function value on the grid, (M+1)*(N+1) matrix.
% x,y:      The (x,y) coordination of the function value in V, matrix with
%           the same size of V. Mainly for the usage of function 'surf.m'.
% Std & Rule:   Please check the function 'quasi_newton.m'.

% Date:     Apr 15th, 2018
% Creator:  BroC


m = 32; n = 32;
type = 6;
u = u_bound(m,n,type);

% set the boundaries
V = zeros(m + 1, n + 1);
V(1,1:n) = u.up;
V(m+1,2:n+1) = u.down;
V(2:m+1,1) = u.left;
V(1:m,n+1) = u.right;

% initial value of v: (1, 1, ..., 1)'
V(2:m,2:n) = ones(m - 1,n - 1);

% coordination
y = ((0:m)'./m) * ones(1,n+1);
x = ones(m+1,1) * ((0:n)./ n);

h = [1 / m, 1 / n];
Std = [1e-6,1e-6,500,inf];
Rule_els.opt = [0, 100, 1, 0.01, 2, 1e-8];

% initialization graph
surf(x,y,V);
v = reshape(V(2:m,2:n), (m - 1) * (n - 1), 1);

% BFGS-H
[vStar, info, perf] = quasi_newton(@fg_vec2surf, [2 1], ...
    v,Std,Rule_els,u,h);

% minimal surface graph
V(2:m,2:n) = reshape(vStar, m - 1, n - 1);
surf(x,y,V);