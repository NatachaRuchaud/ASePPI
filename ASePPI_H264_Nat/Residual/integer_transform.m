function [W]= integer_transform(X)

% X is a 4x4 block of data
% W is the trasnsformed coefficients

a = 1/4;
b = 1/10;
c = sqrt(1/40);

% E is the scaling factor matrix, to make it similar to DCT ?
% Refer to MPEG4-AVC slides (simplified from the H.264 white paper)

E = [a c a c
     c b c b
     a c a c
     c b c b];

 % C is the core transform matrix
C =  [1 1 1 1
      2 1 -1 -2
      1 -1 -1 1
      1 -2 2 -1];
 
W = (C*X*C');
  %W = (double(C)*double(X)*double(C')); 
% W = (C*X*C').*E;  
% TEST NAT FAUX !!!
%  % C is the core transform matrix
% d = 0.5;
% a=0.5;
% b=0.653;
% c=0.271;
% %C =  [1 1 1 1
% %      1 d -d -1
% %      1 -1 -1 1
% %      d -1 1 -d];
% %Ct =  [1 1 1 d
% %      1 d -1 -1
% %      1 -d -1 1
% %      1 -1 1 -d];
% E = [a*a a*b/2 a*a a*b/2
%      a*b/2 b*b/4 a*b/2 b*b/4
%      a*a a*b/2 a*a a*b/2
%      a*b/2 b*b/4 a*b/2 b*b/4];
% 
% W = (C*X*transpose(C));
% W = (C*X*transpose(C)).*E;


% %%%%INVERSE
%  % Ci is the inverse core transform matrix
% Ci =  [1 1 1 1
%       1 1/2 -1/2 -1
%       1 -1 -1 1
%       1/2 -1 1 -1/2];
% 
%  %Y = transpose(Ci)*W*Ci;
% 
% a=0.5;
% b=0.653;
% c=0.271;
% E = [a*a a*b a*a a*b
%      a*b b*b a*b b*b
%      a*a a*b a*a a*b
%      a*b b*b a*b b*b];
%   Y = transpose(Ci)*(W.*E)*Ci;


%A =[a a a a; b c -c -b; a -a -a a; c -b b c];
%At =[a b a c; a c -a -b; a -c -a b; a -b a -c];
%W = (A*X*transpose(A))
%W(abs(W)<10) = 0;
%W = round(W)
%Inv = (transpose(A)*W*A)

end
