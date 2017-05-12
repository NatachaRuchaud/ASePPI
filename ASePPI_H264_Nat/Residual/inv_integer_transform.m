function [Y] = inv_integer_transform(W)

a = 1/4;
b = 1/10;
c = sqrt(1/40);

% E is the scaling factor matrix
% Refer to MPEG4-AVC slides (simplified from the H.264 white paper)

E = [a c a c
     c b c b
     a c a c
     c b c b];

 % Ci is the inverse core transform matrix
Ci =  [1 1 1 1
      1 1/2 -1/2 -1
      1 -1 -1 1
      1/2 -1 1 -1/2];

 %Y = transpose(Ci)*W*Ci;

%%NAT CODE FAUXXX !!
%a=0.5;
%b=0.653;
%c=0.271;
%E = [a*a a*b a*a a*b
%     a*b b*b a*b b*b
%     a*a a*b a*a a*b
%     a*b b*b a*b b*b];

%A =[a a a a; b c -c -b; a -a -a a; c -b b c];
%%%%%%

 Y = Ci'*W*Ci;
  %Y = transpose(Ci)*(W.*E)*Ci;
  %Y = Y/64;

 
end
