function [Z] = quantization(W,QP)

% q is qbits
q = 15 + floor(QP/6);

% M is the multiplying factor which is found from QP value
% MF is the multiplying factor matrix
% rem(QP,6) alpha   beta    gamma
%           (a)     (b)      (g)
% 0         13107   5243    8066
% 1         11916   4660    7490
% 2         10082   4194    6554
% 3         9362    3647    5825
% 4         8192    3355    5243
% 5         7282    2893    4559

MF =[13107 5243 8066
     11916 4660 7490
     10082 4194 6554
     9362  3647 5825
     8192  3355 5243
     7282  2893 4559];
 
x = rem(QP,6);
 
a = MF(x+1,1);
b = MF(x+1,2);
g = MF(x+1,3);

M = [a g a g
     g b g b
     a g a g
     g b g b];

% a = 0.5;
% b = sqrt(2/5);
% g = a*b/2;
% 
% PF = [a*a g a*a g
%      g b*b/4 g b*b/4
%      a g a g
%      g b*b/4 g b*b/4];

% scaling and quantization 
Z = round(W.*(M/(2^q)));



% %%%INVERSE
% 
% % q is qbits
% q = 15 + floor(QP/6);
% %q = 2^floor(QP/6);
% 
% % The scaling factor matrix V depend on the QP and the position of the
% % coefficient.
% %   delta lambda miu
% SM = [10 16 13
%       11 18 14
%       13 20 16
%       14 23 18
%       16 25 20
%       18 29 23];
%  
%  x = rem(QP,6);
%  
%  % find delta, lambda and miu values
%  d = SM(x+1,1);
%  l = SM(x+1,2);
%  m = SM(x+1,3);
% 
%  V = [d m d m
%       m l m l
%       d m d m
%       m l m l];
%   
%  % find the inverse quantized coefficients
%   Wi = Z.*V;
%   %Wi = Wi*q;
%   Wi = bitshift(Wi,q-15,'int64')

 
end
