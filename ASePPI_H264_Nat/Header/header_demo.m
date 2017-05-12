function [bits] = header(h,w,QP, IP, Frame_start,Frame_end, RoI)

bits = '';
bits = [bits dec2bin(h,16) dec2bin(w,16) dec2bin(QP,16) dec2bin(IP,16) dec2bin(Frame_start,16) dec2bin(Frame_end,16)];
s_RoI = size(RoI);
for i=1:s_RoI(1)
 bits = [bits dec2bin(RoI(i, 1), 16) dec2bin(RoI(i, 2), 16) dec2bin(RoI(i, 3), 16) dec2bin(RoI(i, 4), 16)];
end
