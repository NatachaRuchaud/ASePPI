%---------------------------------------------------------
%% Splitting step: Extracted the scrambled coefficients - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was first created by 
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the reference paper (Section 3.4): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function splitData = splitDct(scrambledR)
%Extract the coefficients according to the zigzag scanning
 Z = zigzag(scrambledR);
%Extract only the AC coefficients which are the original scrambled coefficients (DC+AC)
z = Z(2:end);
%Fill the last coefficient to 0 because we lost this information
z(sqrt(16)*sqrt(16)) = 0;
splitData= z;
splitData = izigzag(splitData, sqrt(16), sqrt(16));

Zeros = zeros(sqrt(16), sqrt(16));
Zeros(1:sqrt(16), 1:sqrt(16)) = splitData;
splitData = Zeros;
end

