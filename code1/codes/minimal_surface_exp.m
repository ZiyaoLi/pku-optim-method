% MINIMAL_SURFACE_EXP.M solves mininal surface with boundary constraints
% problem with differential method. Detailed description is in "p153",
% Problem 4.3. In this script, the convex set D is limited to a
% square (rectangular) grid, i.e. (0,1)*(0,1) \in R^2. Numerical
% Experiments are conducted in this script.
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
for type=0:6
    u = u_bound(m,n,type);
    
    % initial value of v: (1, 1, ..., 1)'
    v = ones((m - 1) * (n - 1), 1);
    h = [1 / m, 1 / n];
    Std = [1e-6,1e-6,500,inf];
    Rule_els.opt = [0, 100, 1, 0.01, 2, 1e-8];
        
    % Damped Newton:
    [~, info, ~] = damp_newton(@fg_vec2surf, @eval_diff_G, ...
        v,Std,Rule_els,2,u,h,@fg_vec2surf,1e-8,u,h);
    Record.iter(type + 1, 1) = info(2);
    Record.feva(type + 1, 1) = info(3);
    Record.stat(type + 1, 1) = info(1);
    
    % LM-Adjusted Newton:
    [~, info, ~] = lm_newton(@fg_vec2surf, @eval_diff_G, ...
        v,Std,Rule_els,2,u,h,@fg_vec2surf,1e-8,u,h);
    Record.iter(type + 1, 2) = info(2);
    Record.feva(type + 1, 2) = info(3);
    Record.stat(type + 1, 2) = info(1);
    
    % SR1-B:
    [~, info, ~] = quasi_newton(@fg_vec2surf,[0 1],...
        v,Std,Rule_els,u,h);
    Record.iter(type + 1, 3) = info(2);
    Record.feva(type + 1, 3) = info(3);
    Record.stat(type + 1, 3) = info(1);
    
    % SR1-H:
    [~, info, ~] = quasi_newton(@fg_vec2surf,[0 0],...
        v,Std,Rule_els,u,h);
    Record.iter(type + 1, 4) = info(2);
    Record.feva(type + 1, 4) = info(3);
    Record.stat(type + 1, 4) = info(1);
    
    % DFG-B:
    [~, info, ~] = quasi_newton(@fg_vec2surf,[1 1],...
        v,Std,Rule_els,u,h);
    Record.iter(type + 1, 5) = info(2);
    Record.feva(type + 1, 5) = info(3);
    Record.stat(type + 1, 5) = info(1);
    
    % DFG-H:
    [~, info, ~] = quasi_newton(@fg_vec2surf,[1 0],...
        v,Std,Rule_els,u,h);
    Record.iter(type + 1, 6) = info(2);
    Record.feva(type + 1, 6) = info(3);
    Record.stat(type + 1, 6) = info(1);
    
    % BFGS-B:
    [~, info, ~] = quasi_newton(@fg_vec2surf,[2 1],...
        v,Std,Rule_els,u,h);
    Record.iter(type + 1, 7) = info(2);
    Record.feva(type + 1, 7) = info(3);
    Record.stat(type + 1, 7) = info(1);
    
    % BFGS-H:
    [~, info, ~] = quasi_newton(@fg_vec2surf,[2 0],...
        v,Std,Rule_els,u,h);
    Record.iter(type + 1, 8) = info(2);
    Record.feva(type + 1, 8) = info(3);
    Record.stat(type + 1, 8) = info(1);
end