function [ f, g ] = fg_test( x, calc_g, p )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

f = x * x + p;
g = false;
if calc_g
    g = 2 * x;
end
end

