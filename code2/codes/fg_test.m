function [ f, g ] = fg_test( x, calc_g, p )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

f = x * x + p;
g = false;
if calc_g
    g = 2 * x;
end
end

