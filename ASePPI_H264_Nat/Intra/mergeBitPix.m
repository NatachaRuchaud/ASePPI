%---------------------------------------------------------
%% Merging step: Hidding the scrambled coefficients inside the AC ones while inserting the DC chosen -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was first created by
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the reference paper (Section 3.2.3): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mergeData = mergeBitPix(dct_Edg_quant, scrambledR)
%Extract the coefficients according to the zigzag scanning
zs=zigzag(scrambledR);
%%Insert the DC chosen into the first coefficient
Z = [dct_Edg_quant(1, 1), zs];
z = zeros(1, 64);
z(1:length(Z)) = Z;
mergeData =izigzag(z, sqrt(16), sqrt(16));
end

