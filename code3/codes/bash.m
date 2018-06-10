% BASH.M solves six least square problems mentioned in //Testing
% Unconstrained Optimization Software//, namely Meyer function, namely
% Meyer Function, Osborne Exponential I Function, Osborne Exponential II
% Function, Jennrich & Sampson Funtion, Freudenstein and Roth Function and 
% Davidon Function function. Six different methods are examined. Generally 
% Gauss-Newton (GN) method have higher convergence speed, while the exact 
% line search used in it adds greatly to its function-evaluation times. 
% Levenberg-Marquardt-Fletcher (LMF) method and Levenberg-Marquardt-Nielson 
% (LMN) method show similar performances. Dogleg and Twofold Dogleg methods 
% do not perform well on the problem. Gauss-Newton method with the BFGS 
% adjustment for Large-Residual problem is the most robust method with more 
% iterations to converge in different problems.

% Arguments
% frjs & x0s:   Function handles and initial points of corresponding
%               problems.
% Std & Rule:   Please check the relavent optimization functions.

% Date:     Jun 11th, 2018
% Creator:  BroC

Std = [1e-8,1e-8,1000,Inf];
Rule.opt = [0, 100, 1, 0.1, 2, 1e-8];
iterations = zeros(3, 6, 6);
frjs = {@frj_meyer, @frj_osborne1, @frj_osborne2, ...
        @frj_js,    @frj_fr,       @frj_davidon};
x0s = {[.02, 4000, 250]', ...
       [.5, 1.5, -1, .01, .02]', ...
       [1.3, .65, .65, .7, .6, 3, 5, 7, 2, 4.5, 5.5]', ...
       [.3, .4]', ...
       [.5, -2]', ...
       [25, 5, -5, 1]'};
rs = cell(6);

for i=1:6
    fprintf('Problem %d - Gauss Newton\n', i);
    [~, info, perf] = gauss_newton(frjs{i}, x0s{i}, Std, Rule);
    iterations(1:3, i, 1) = info;
    rs{1, i} = perf;

    fprintf('Problem %d - LMF\n', i);
    [~, info, perf] = lm(frjs{i}, 1, x0s{i}, Std);
    iterations(1:3, i, 2) = info;
    rs{2, i} = perf;

    fprintf('Problem %d - LMN\n', i);
    [~, info, perf] = lm(frjs{i}, 2, x0s{i}, Std);
    iterations(1:3, i, 3) = info;
    rs{3, i} = perf;
    
    fprintf('Problem %d - Dogleg\n', i);
    [~, info, perf] = dogleg(frjs{i}, 1, x0s{i}, Std, 2);
    iterations(1:3, i, 4) = info;
    rs{4, i} = perf;
    
    fprintf('Problem %d - Twofold Dogleg\n', i);
    [~, info, perf] = dogleg(frjs{i}, 2, x0s{i}, Std, 2);
    iterations(1:3, i, 5) = info;
    rs{5, i} = perf;

    fprintf('Problem %d - Gauss Newton + BFGS \n', i);
    [~, info, perf] = gauss_bfgs(frjs{i}, x0s{i}, Std, Rule);
    iterations(1:3, i, 6) = info;
    rs{6, i} = perf;
end

