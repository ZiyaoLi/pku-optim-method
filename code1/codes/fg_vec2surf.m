function [ f, g ] = fg_vec2surf( v,u,h )
% FG_VEC2SURF.M function calculates the area of the surface based on grid
% values provided in v & u, and the grid length h.
%
% Input:
% v:        The vectorized matrix of the values on the grid V.
%           (M-1)*(N-1) vector, where M,N are the number of grids of V.
% u:        The boundary constraints generated from function u_bound.m,
%           struct.
% h:        The different grid length on the grids, 2-vector.
%           h(1) = hx, h(2) = hy.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = fg_vec2surf(v, u, h)

% Date:     Apr 14th, 2018
% Creator:  BroC

% pre-calculation to save time
m = length(u.left);
n = length(u.up);
txy = h(1) * h(1) * h(2) * h(2);
tx = h(1) * h(1); ty = h(2) * h(2);

% forming the actual grids V
V = zeros(m + 1, n + 1);
V(2:m,2:n) = reshape(v, m - 1, n - 1);
V(1,1:n) = u.up;
V(m+1,2:n+1) = u.down;
V(2:m+1,1) = u.left;
V(1:m,n+1) = u.right;

% initialization
f = 0; g = zeros(m + 1, n + 1);

% calculation
for i = 1:m
    for j = 1:n
        dd = V(i+1, j) - V(i,j);
        dr = V(i, j+1) - V(i,j);
        fL = sqrt(txy + ty * dd^2 + tx * dr^2) * .5;
        f = f + fL;
        g(i+1, j) = g(i+1, j) + .25 * ty * dd / fL;
        g(i, j+1) = g(i, j+1) + .25 * tx * dr / fL;
        g(i, j) = g(i, j) - .25 * (ty * dd + tx * dr) / fL;
    end
end
for i=2:m+1
    for j=2:n+1
        du = V(i-1, j) - V(i,j);
        dl = V(i, j-1) - V(i,j);
        fU = sqrt(txy + ty * du^2 + tx * dl^2) * .5;
        f = f + fU;
        g(i-1, j) = g(i-1, j) + .25 * ty * du / fU;
        g(i, j-1) = g(i, j-1) + .25 * tx * dl / fU;
        g(i, j) = g(i, j) - .25 * (ty * du + tx * dl) / fU;
    end
end
g = reshape(g(2:m,2:n), (m - 1) * (n - 1), 1);
end

