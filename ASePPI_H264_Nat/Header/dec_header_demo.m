function [h,w,QP,IP, Frame_start,Frame_end,m, RoI] = dec_header(bits)

m = 1;

h = bin2dec(bits(m:m+15));
m = m + 16;
disp('Height=')
disp(h)

w = bin2dec(bits(m:m+15));
m = m + 16;
disp('Width=')
disp(w)

QP = bin2dec(bits(m:m+15));
m = m + 16;
disp('QP=')
disp(QP)

IP = bin2dec(bits(m:m+15));
m = m + 16;
disp('IP=')
disp(IP)

Frame_start = bin2dec(bits(m:m+15));
m = m + 16;
disp('Frame_start=')
disp(Frame_start)

Frame_end = bin2dec(bits(m:m+15));
m = m + 16;
disp('Frame_end=')
disp(Frame_end)

for i = 1:ceil(Frame_end/IP)
RoI(i, 1) = bin2dec(bits(m:m+15));
m = m + 16;
RoI(i, 2) = bin2dec(bits(m:m+15));
m = m + 16;
RoI(i, 3) = bin2dec(bits(m:m+15));
m = m + 16;
RoI(i, 4) = bin2dec(bits(m:m+15));
m = m + 16;
end
disp('RoI=')
disp(RoI)
