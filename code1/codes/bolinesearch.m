function [alpha, info, perf] = ...
    bolinesearch(fg_handle, x0, dir, Rule, varargin)
% BOLINESEARCH Find the answer to problem P: a = argmin f(x+a*d)
%
% This function will use exact or inexact line search to solve the problem.
% When doing inexact line search, you can choose Armijo-Goldstein, Wolfe,
% strong Wolfe or you own criterion.
% Denote g as the gradient of f. 
%
% Input 
% fg_handle:	f & g in P, function_handle
%               The function should be declared as a M-function
%                   [F g] = f(x, p1,p2,...)
%               F is a scalar, along with g an n-vector.
% x0:           x0 in P, n-vector
% dir:          d in P, n-vector
% Rule:     option & method & criterion to solve P, struct
%       Rule.crtr:	criterion, function_handle
%               The function should be declared as a M-function
%                   Judge = criterion(Step, StepSize, Data0, Data, flag).
%               Judge is a logical number with 1 perfect, and Step is just
%               as mentioned above StepSize below, and Data0 & Data are
%               structs including possible fields F g, flag is a array
%               including several parameters in the criterion.
%                - options: boarmgld, bowlf, bostwlf
%       Rule.mthd:	method to get new point, function_handle
%               The function shoule be declared as a M-function
%                   [NewSize New] = method(ObjFun, Point, Step, StepSize,  
%                                          Data, p1,p2,...).
%               Notice here StepSize & NewSize can be a scalar array(to
%               make coding easy), but only the 1st element is just needed
%               actually. And in that situation the Data & New will be
%               struct arrays.
%                - options: bointrplt22, bointrplt33
%       Rule.opt:       options of iteration, scalar array
%           opt(1):     0 - exact line search, use .618 method
%                       1 - inexact, need crtr & mthd
%              (2): maximum of iterations
%              (3): upper bound of a
%            (4:5): criterion flag, but also can be .618 method flag which
%                   have one element e with default 1e-3
%            (3:6): alpha0, gamma0, t, eps for ELS (0.618 method)
%       default:        bostwlf - bointrplt33 - [1 10 10 0.95 0.05]
%                       or [0 10 25 1e-4] (if opt(1) = 0)
%
% Output
% StepSize: a in P just the answer, scalar
% info:     status to execute this function, integral array
%       info(1):    exit code
%               0 - Successful call
%               1 - Reach the maximum of iterations
%              -1 - a is not legal
%           (2):    number of iterations
%           (3):    number of evaluating ObjFun
% perf:     other useful data, struct
%       perf.x:     new point after iteration
%       perf.F:     function value at new point
%       perf.g:     gradient at new point
%
% Call
% [StepSize err] = bolinesearch(ObjFun, Point, Step)
% [StepSize err] = bolinesearch(ObjFun, Point, Step, Rule)
% [StepSize err] = bolinesearch(ObjFun, Point, Step, Rule, p1,p2,...)
% [StepSize err perf] = bolinesearch(......)

% Version:  2009.04.25
% Create:   2009.04.23
% Coder:    Xin Liang
% Update:   Apr 14th, 2018
% Updater:  broC


narginchk(3, inf);
nargoutchk(0, 3);
if nargin == 3 || isempty(Rule)
    Rule.opt = [0 1 0.1 2 1e-8];
end

info = [0 0 0];
[Point0.F, Point0.g] = feval(fg_handle, x0, varargin{:});

if ~Rule.opt(1) % Execute ELS process
    
    k = ( sqrt(5) - 1 ) / 2; % 0.618
    
    % search for an eligible region [a,b]
    
    while true
        a = Rule.opt(3);
        gm = Rule.opt(4);
        t = Rule.opt(5);
        eps = Rule.opt(6);
        fp = feval(fg_handle, x0 + a .* dir, varargin{:});
        info(3) = info(3) + 1;
        while fp == Inf
            % indicating overflow; cutting down the alpha
            a = a / t;
            gm = gm / t;
            fp = feval(fg_handle, x0 + a .* dir, varargin{:});
            info(3) = info(3) + 1;
        end
        flag_first = true;
        while true
            ap = a;
            a = a + gm;
            f = feval(fg_handle, x0 + a .* dir, varargin{:});
            info(3) = info(3) + 1;
            if f >= fp || a <= 0
                if a <= 0
                    a = 0;
                end
                if flag_first
                    gm = -gm;
                    alpha = a;
                    a = ap;
                    flag_first = false;
                    continue
                else
                    L = min(alpha, a);
                    U = max(alpha, a);
                    break
                end
            else
                gm = t * gm;
                alpha = ap;
                flag_first = false;
            end
        end

        lb=L; ub=U;
        fl = false;
        fr = false;
        while ub - lb >= eps
            if ~fl
                al = lb + (1-k) * (ub - lb); % 0.382 point calculated
                Fl = feval(fg_handle, x0 + al * dir, varargin{:});
                info(3) = info(3) + 1;
            end
            if ~fr
                ar = lb + k * (ub - lb); % 0.618 point calculated
                Fr = feval(fg_handle, x0 + ar * dir, varargin{:});
                info(3) = info(3) + 1;
            end
            if Fl < Fr
                ub = ar;
                ar = al; Fr = Fl;
                fr = true; fl = false;
            else
                lb = al;
                al = ar; Fl = Fr;
                fr = false; fl = true;
            end
            info(2) = info(2) + 1;
            if info(2) >= Rule.opt(2)
                info(1) = 1;
                break;
            end
        end
        alpha = (lb + ub) / 2;
        perf.x = x0 + alpha .* dir;
        [perf.F, perf.g] = feval(fg_handle, perf.x, varargin{:});
        if Point0.F > perf.F
            break;
        elseif alpha < 1e-8
            % which means no step length is derived.
            fprintf(['Oops, No step length found in the ELS searching.\n'...
                    'Consider changing a searching method.\n'])
            N_sample = 50;
            x_plot = (0:N_sample) ./ N_sample .* U;
            y_plot = zeros(1, N_sample + 1);
            for i=1:(N_sample + 1)
                y_plot(i) = feval(fg_handle, x0 + x_plot(i) .* dir, varargin{:});
            end
            plot(x_plot,y_plot)
            perf.x = x0;
            perf.F = Point0.F;
            perf.g = Point0.g;
            alpha = 0;
            info(1) = -1;
            break;
        else
            % which means the ELS result is even larger than f(x0).
            Rule.opt(3) = Rule.opt(3) / t;
            Rule.opt(4) = Rule.opt(4) / t;
        end
    end
    
else % Execute ILS process
   [Point1.F, Point1.g] = feval(fg_handle, x0 + Rule.opt(3) .* dir,...
       varargin{:});
   info(3) = info(3) + 2;
   Judge = feval(Rule.crtr, dir, Rule.opt(3), Point0, Point1, Rule.opt(4:5));
   while ~Judge
       [alpha, Point1] = feval(Rule.mthd, ...
           fg_handle, Point1, dir, alpha, Point1, varargin{:});
       info(2) = info(2) + 1;
       info(3) = info(3) + 1;
       if ~isreal(alpha(1))
           % cannot derive a real alpha in the interpolation process.
           info(1) = -1;
           break;
       end
       if alpha(1) < 0
           % cannot derive a positive alpha in the interpolation process.
           info(1) = -1;
           alpha(1) = 0;
       end
       Judge = feval(Rule.crtr, ...
           dir, alpha(1), Point0, Point1(1), Rule.opt(4:5));
       if info(2) >= Rule.opt(2)
           info(1) = 1;
           break;
       end
   end
   alpha = alpha(1);
   perf.x = Point1 + alpha*dir;
   perf.F = Point1(1).F;
   perf.g = Point1(1).g;
end
