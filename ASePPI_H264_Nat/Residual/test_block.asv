
clear all;
clc;

load table;
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros

Seq= [200 -1 -3 1 0 0 0 0
     -7 -4 -5 2 0 0 0 0
     -1 -6 2 4 0 0 0 0
     1 -1 20 1 0 0 0 0
     -2 4 0 -1 0 0 1 0
     3 0 0 0 0 0 0 0
     -3 0 0 0 1 0 0 1
     0 0 0 0 -1 0 0 0];
 
 QP = 10;
 bits = '';
 [h,w]=size(Seq);
 
for i=1:4:h
    for j=1:4:w
        X = Seq(i:i+3,j:j+3);
%         W = integer_transform(X);
%         Z = quantization(W,QP);
        [bits_new] = enc_cavlc(X, 0, 0);
        bits = [bits bits_new];
    end
end

k=1;
while (k<=length(bits))
    for i=1:4:h
        for j=1:4:w
            [data,m] = dec_cavlc(bits(k:length(bits)),0,0);
            Seq_r(i:i+3,j:j+3) = data;
            k = k + m - 1;
        end
    end
end

diff = Seq - Seq_r
%  Wi = inv_quantization(Z,QP);
% 
%  Y = inv_integer_transform(Wi);
% 
%  %  post scaling - very important 
%  Xi = round(Y/64);
 
 