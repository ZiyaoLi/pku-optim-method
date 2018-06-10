% BASH.M solves three optimization problems mentioned in Extended CUTEr, 
% namely generalized BARD function, generalized CUBE function and s303-305 
% function. BB and other 5 conjugate gradient methods are attempted, namely
% FR (Fletcher-Reeves), PRP (Polak-Ribiere-Polyak), PRP+ (Polak-Ribiere-
% Polyak Plus), CD (Conjugate Descend) and DY (Dai-Yuan). Both iteration 
% and function evaluation standards are examined and saved to compare 
% performance.
% 
% Arguments
% list_n:   The candidate scale of Problem Box Three-Dimensional
%           Function (Penalty II Function).
% x0_box (x0_pen):  The initial values of the two problem w.r.t. problem
%           scale.
% Std & Rule:   Please check the relavent optimization functions.

% Date:     Apr 13th, 2018
% Creator:  BroC

list_n = [1e2, 1e3, 1e4, 1e5];

Std = [1e-2,1e-2,1000,Inf];
Rule_els.opt = [0, 100, 2, 0.1, 2, 1e-8];
Rule_els_s303.opt = [0, 300, 1e-8, 1e-9, 10, 1e-40];
Rule_els_bard.opt = [0, Inf, 1e-2, 5e-3, 2, 1e-12];
Args_BB = [1 5 .8 1 .5 1e-8];
Args_s303 = [1 5 .8 1 .5 1e-40];
y = [14 18 22 25 29 32 35 39 37 58 73 96 134 210 439] ./ 100;
iterations = zeros(3, 6, 3, 3);

for i=1:3
    x0 = ones(list_n(i), 1);
    x0_s303 = x0 ./ 10;
    x0_cube = x0 .* 0.9;
    for j=1:5
        fprintf('Gen-BARD  n=%3d  method=%d\n', list_n(i), j);
        [~, info, ~] = conjugate_gradient(@fg_genbard, j, ...
            x0, Std, Rule_els_bard, y);
        iterations(1:3, j, 1, i) = info;
        
        fprintf('Gen-CUBE  n=%3d  method=%d\n', list_n(i), j);
        [~, info, ~] = conjugate_gradient(@fg_gencube, j, ...
            x0_cube, Std, Rule_els);
        iterations(1:3, j, 2, i) = info;
                
        fprintf('s303-305  n=%3d  method=%d\n', list_n(i), j);
        [~, info, ~] = conjugate_gradient(@fg_s303, j, ...
            x0_s303, Std, Rule_els_s303);
        iterations(1:3, j, 3, i) = info;
    end
    
    fprintf('Gen-BARD  n=%3d  method=BB\n', list_n(i));
    [~, info, ~] = global_bb(@fg_genbard, Args_BB, ...
            x0, Std, Rule_els, y);    
    iterations(1:3, 6, 1, i) = info;

    fprintf('Gen-CUBE  n=%3d  method=BB\n', list_n(i));
    [~, info, ~] = global_bb(@fg_gencube, Args_BB, ...
            x0_cube, Std, Rule_els);
    iterations(1:3, 6, 2, i) = info;
    
    fprintf('s303-305  n=%3d  method=BB\n', list_n(i));
    [~, info, ~] = global_bb(@fg_s303, Args_s303, ...
            x0_s303, Std, Rule_els);    
    iterations(1:3, 6, 3, i) = info;
    
end

