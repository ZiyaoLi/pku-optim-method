function opt = bodfltfunc(opt, default)
%BODLTCHK Fill default function_handle to opt when blank
%
% Call
% optcrct = bodfltfunc(opt, default)

% Version:  2009.04.25
% Create:   2009.04.24
% Coder:    Xin Liang


narginchk(2, 2);
nargoutchk(0, 1);

if ~isa(opt, 'function_handle')
    opt = default;
end
