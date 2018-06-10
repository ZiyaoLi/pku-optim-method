% BASH.M solves Box three-dimensional function and Penalty II function
% with different kinds of optimization algorithms, namely damped newton 
% method, fixed newton method (LM method) and quasi-newton methods. Both 
% iteration and function evaluation standards are examined and saved to 
% compare performance.
% 
% Arguments
% list_m (list_n):  The candidate scale of Problem Box Three-Dimensional
%           Function (Penalty II Function).
% x0_box (x0_pen):  The initial values of the two problem w.r.t. problem
%           scale.
% Std & Rule:   Please check the relavent optimization functions.

% Date:     Apr 13th, 2018
% Creator:  BroC


list_m = [3,5,10,15,20];
list_n = [2,4,6,8,10];
x0_box = [0,10,20]';
x0_pen = cell(5,1);
for t=1:5
    x0_pen{t} = 0.5.*ones(list_n(t),1);
end

Std = [1e-8,1e-8,2000,inf];
Rule_els.opt = [0, 100, 2, 0.1, 2, 1e-8];
Rule_ils.opt = [1, 100, 100, 0.2, 0.1];
Rule_ils.crtr = @bowlf;
Rule_ils.mthd = @bointrplt22;


% Box 3D problem
for i=1:5
    % Damped Newton:
    [~, info, ~] = damp_newton(@fg_box3d, @eval_diff_G, ...
        x0_box,Std,Rule_els,1,list_m(i),@fg_box3d,1e-8,list_m(i));
    Record_box.iter(i,1) = info(2);
    Record_box.feva(i,1) = info(3);
    Record_box.stat(i,1) = info(1);
    
    % LM-Adjusted Newton:
    Std_special = [1e-4,1e-4,10000,inf];
    [~, info, ~] = lm_newton(@fg_box3d, @eval_diff_G, ...
        x0_box,Std,Rule_els,1,list_m(i),@fg_box3d,1e-8,list_m(i));
    Record_box.iter(i,2) = info(2);
    Record_box.feva(i,2) = info(3);
    Record_box.stat(i,2) = info(1);
    % result: iter = |...|8;381(1e-4)|...|...|...|
    % Cannot derive step lengths when precision is higher
    % (Maximum Iter Reached).
    
    % SR1-B:
    [~, info, ~] = quasi_newton(@fg_box3d, [0 1], ...
        x0_box,Std,Rule_els,list_m(i));
    Record_box.iter(i,3) = info(2);
    Record_box.feva(i,3) = info(3);
    Record_box.stat(i,3) = info(1);
    % result: iter = |...|3;153(1e-1)|3;163(1e0)|3;158(1e0)|...|
    % Cannot derive step lengths when precision is higher
    % (Wrong Directions).
    
    % SR1-H:
    [~, info, ~] = quasi_newton(@fg_box3d, [0 0], ...
        x0_box,Std,Rule_els,list_m(i));
    Record_box.iter(i,4) = info(2);
    Record_box.feva(i,4) = info(3);
    Record_box.stat(i,4) = info(1);
    % result: iter = |...|3;153(1e-1)|3;163(1e0)|3;158(1e0)|...|
    % Cannot derive step lengths when precision is higher
    % (Wrong Directions).
    
    % DFG-B:
    [~, info, ~] = quasi_newton(@fg_box3d, [1 1], ...
        x0_box,Std,Rule_els,list_m(i));
    Record_box.iter(i,5) = info(2);
    Record_box.feva(i,5) = info(3);
    Record_box.stat(i,5) = info(1);
    
    % DFG-H:
    [~, info, ~] = quasi_newton(@fg_box3d, [1 0], ...
        x0_box,Std,Rule_els,list_m(i));
    Record_box.iter(i,6) = info(2);
    Record_box.feva(i,6) = info(3);
    Record_box.stat(i,6) = info(1);
    
    % BFGS-B:
    [~, info, ~] = quasi_newton(@fg_box3d, [2 1], ...
        x0_box,Std,Rule_els,list_m(i));
    Record_box.iter(i,7) = info(2);
    Record_box.feva(i,7) = info(3);
    Record_box.stat(i,7) = info(1);
    
    % BFGS-H:
    [~, info, ~] = quasi_newton(@fg_box3d, [2 0], ...
        x0_box,Std,Rule_els,list_m(i));
    Record_box.iter(i,8) = info(2);
    Record_box.feva(i,8) = info(3);
    Record_box.stat(i,8) = info(1);
end


% Penalty II problem
for i=1:5
    % Damped Newton:
    [~, info, ~] = damp_newton(@fg_penalty2, @eval_diff_G, ...
        x0_pen{i},Std,Rule_els,1,list_n(i),@fg_penalty2,1e-8,list_n(i));
    Record_pen.iter(i,1) = info(2);
    Record_pen.feva(i,1) = info(3);
    Record_pen.stat(i,1) = info(1);
    % result: iter = | - |4;192(1e-4)|3;137(1e-5)|...|4;190(1e-3)|
    % Cannot derive step lengths when precision is higher
    % (Wrong Directions).
    
    % LM-Adjusted Newton:
    Std_special = [1e-4,1e-4,10000,inf];
    [~, info, ~] = lm_newton(@fg_penalty2, @eval_diff_G, ...
        x0_pen{i},Std,Rule_els,1,list_n(i),@fg_penalty2,1e-8,list_n(i));
    Record_pen.iter(i,2) = info(2);
    Record_pen.feva(i,2) = info(3);
    Record_pen.stat(i,2) = info(1);
    % result: iter = |...|8;381(1e-4)|...|...|...|
    % Cannot derive step lengths when precision is higher
    % (Maximum Iter Reached).

    % SR1-B:
    [~, info, ~] = quasi_newton(@fg_penalty2, [0 1], ...
        x0_pen{i},Std,Rule_els,list_n(i));
    Record_pen.iter(i,3) = info(2);
    Record_pen.feva(i,3) = info(3);
    Record_pen.stat(i,3) = info(1);
    % result: iter = |1;49(1e0)| - |7;339(1e-5)|3;158(1e0)|2;97(1e-1)|
    % Cannot derive step lengths when precision is higher
    % (Wrong Directions).

    % SR1-H:
    [~, info, ~] = quasi_newton(@fg_penalty2, [0 0], ...
        x0_pen{i},Std,Rule_els,list_n(i));
    Record_pen.iter(i,4) = info(2);
    Record_pen.feva(i,4) = info(3);
    Record_pen.stat(i,4) = info(1);
    % result: iter = |1;49(1e0)| - |7;339(1e-5)|3;158(1e0)|2;97(1e-1)|
    % Cannot derive step lengths when precision is higher
    % (Wrong Directions).

    % DFP-B:
    [~, info, ~] = quasi_newton(@fg_penalty2, [1 1], ...
        x0_pen{i},Std,Rule_els,list_n(i));
    Record_pen.iter(i,5) = info(2);
    Record_pen.feva(i,5) = info(3);
    Record_pen.stat(i,5) = info(1);
    
    % DFP-H:
    [~, info, ~] = quasi_newton(@fg_penalty2, [1 0], ...
        x0_pen{i},Std,Rule_els,list_n(i));
    Record_pen.iter(i,6) = info(2);
    Record_pen.feva(i,6) = info(3);
    Record_pen.stat(i,6) = info(1);
    
    % BFGS-B:
    [~, info, ~] = quasi_newton(@fg_penalty2, [2 1], ...
        x0_pen{i},Std,Rule_els,list_n(i));
    Record_pen.iter(i,7) = info(2);
    Record_pen.feva(i,7) = info(3);
    Record_pen.stat(i,7) = info(1);
    
    % BFGS-H:
    [~, info, ~] = quasi_newton(@fg_penalty2, [2 0], ...
        x0_pen{i},Std,Rule_els,list_n(i));
    Record_pen.iter(i,8) = info(2);
    Record_pen.feva(i,8) = info(3);
    Record_pen.stat(i,8) = info(1);
end

for i=2:2
    % LM-Adjusted Newton:
    Std_special = [1e-4,1e-4,10000,inf];
    [~, info, ~] = lm_newton(@fg_box3d, @eval_diff_G, ...
        x0_box,Std,Rule_els,1,list_m(i),@fg_box3d,1e-8,list_m(i));
    Record_box.iter(i,2) = info(2);
    Record_box.feva(i,2) = info(3);
    Record_box.stat(i,2) = info(1);
    % result: iter = |...|8;381(1e-4)|...|...|...|
    % Cannot derive step lengths when precision is higher
    % (Maximum Iter Reached).
end
