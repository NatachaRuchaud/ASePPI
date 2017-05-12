%---------------------------------------------------------
%% Scrambling step -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was first created by
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the reference paper (Section 3.3): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs are:
% dct_imgCache_quant - array containing the residual coefficients after the transformation and the quantization
% seed - the seed to generate random numbers
%
% outputs are:
% scramble - the scrambled coefficients

function [scramble] = scrambledXoR_Inter_Forward( dct_imgCache_quant, seed)
%Extract the coefficients according to the zigzag scanning
zigzagScanning = zigzag(dct_imgCache_quant);
zigzagScanningKept = zigzagScanning(1:16);
DC = de2bi(abs(zigzagScanningKept(1)), 16);

rng(seed)
%Extract the DC value and convert to the binary code
dc = abs(zigzagScanningKept(1));
random = round(rand(1)*127);

%%%%Calculate the number of combinaison
%Permutation !!!
npermute1=0;
ind_nn_zeros = find(zigzagScanningKept~=0);
if(length(ind_nn_zeros)~=0)
    zigzagNNzeros = zigzagScanningKept(1:ind_nn_zeros(length(ind_nn_zeros)));
    npermute1 = length(zigzagNNzeros)-1;
    p = randperm(npermute1);
end
%Flipping !!!
ac_nn_zeros = find(zigzagScanning~=0);
if(length(ac_nn_zeros)>=1)
    npermute2 = length(ac_nn_zeros)-1;
else
    npermute2 = length(ac_nn_zeros);
end

%%%%Choosing the method which has the higher number of combinaison
if(gamma(npermute1+1)>=2^npermute2)
    npermute = npermute1;
    method = 1;
else
    npermute =npermute2;
    method = 2;
end

%%Begin of the DC Encryption
if(dc~=0)
    if(dc<16)
        cryptDC = bi2de(xor(DC, de2bi(mod(random, 16), 16)));
        decryptDC = bi2de(xor(de2bi(cryptDC,16), de2bi(mod(random, 16), 16)));
        if(decryptDC~=dc)
            dc
            mod(random, 16)
            cryptDC
            decryptDC
        end
    else
        module = floor(log2(dc));
        cryptDC = bi2de(xor(DC, de2bi(mod(random, 2^(module)), 16)));
    end
    if(cryptDC ==0)
        cryptDC=zigzagScanningKept(1);
    else
        if(zigzagScanningKept(1)<0)
            cryptDC = -cryptDC;
        end
    end
else
    cryptDC = zigzagScanningKept(1);
end
%%End of DC Encryption


%Inverse coefficients signs randomly
if(method==2)
    Rr = round(rand(4, 4));
    [indx_nonzeros, indy_nonzeros] = find(Rr);
    invSigneRandom = dct_imgCache_quant;
    invSigneRandom(indx_nonzeros(:), indy_nonzeros(:)) = -dct_imgCache_quant(indx_nonzeros(:), indy_nonzeros(:));
    invSigneRandom(1) = cryptDC;
    scramble = invSigneRandom;
end

%Permute randomly the coefficients
if(method==1)
    zigzagNNzeros(1) = cryptDC;
    if(npermute~=0);
        permuteZigzag = zigzagNNzeros(p);
        permuteZigzag(npermute+1) = zigzagNNzeros(npermute+1);
        final(1:npermute+1) =  permuteZigzag;
        final(npermute+2:16) = 0;
    else
        final(1)= cryptDC;
        final(2:16) = 0;
    end
    scramble = izigzag(final, sqrt(16), sqrt(16));
end
end


