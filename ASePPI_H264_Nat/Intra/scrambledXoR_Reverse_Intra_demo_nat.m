%---------------------------------------------------------
%% DeScrambling step - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was first created by 
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the reference paper (Section 3.4): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs are:
% dct_imgCache_quant - array containing the residual coefficients after the transformation and the quantization
% DC_Vizu - the DC chosen for the final vizualization
%
% outputs are:
% descramble - the descrambled coefficients

function [descramble] = scrambledXoR_Reverse_Intra_demo_nat( dct_imgCache_quant, DC_Vizu)

%Extract the coefficients according to the zigzag scanning
            zigzagScanning = zigzag(dct_imgCache_quant);
            zigzagScanningKept = zigzagScanning(1:15);
            random = round(rand(1)*127);

%%%%Calculate the number of combinaison
%Permutation !!!
npermute1=0;
ind_nn_zeros = find(zigzagScanningKept~=0);
if(length(ind_nn_zeros)~=0)
    zigzagNNzeros = zigzagScanningKept(1:ind_nn_zeros(length(ind_nn_zeros)));
    npermute1 = length(zigzagNNzeros)-1;
	%If the DC is as the same place
    %npermute1 = npermute1-1;
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
%Inverse AC coefficients signs randomly
if(method==2)
            Rr = round(rand(4, 4));
            [indx_nonzeros, indy_nonzeros] = find(Rr);
            invSigneRandom = dct_imgCache_quant;
            invSigneRandom(indx_nonzeros(:), indy_nonzeros(:)) = -dct_imgCache_quant(indx_nonzeros(:), indy_nonzeros(:));
	    invSigneRandom(1,1) = dct_imgCache_quant(1,1);
            final =invSigneRandom;
end

%Permute randomly the coefficients (DC+AC)
if(method==1)
	if(npermute1~=0)
                        permuteZigzag(p) = zigzagNNzeros(1:npermute);
                        permuteZigzag(npermute+1) = zigzagNNzeros(npermute+1);
			final(1:npermute+1) =  permuteZigzag;
                        final(npermute+2:16) = 0;                  
                else
		    final(1) =zigzagScanning(1);
                    final(2:16)=0;
                end
end
%%Begin of the DC Encryption 
%Extract the DC value and convert to the binary code
                DC = de2bi(abs(final(1)), 16);
		dc = abs(final(1));
                if(dc~=0)
			if(dc<16)
				cryptDC = bi2de(xor(DC, de2bi(mod(random, 16), 16)));
		    else
				module = floor(log2(dc));
				cryptDC = bi2de(xor(DC, de2bi(mod(random, 2^(module)), 16)));
			end
				if(cryptDC ==0)
				    cryptDC=final(1);
				else
				    if(final(1)<0)
					cryptDC = -cryptDC;
				    end
				end
                else
                    cryptDC = final(1);
                end
		final(1) = cryptDC;
		%Decrypt the difference between the DC for the visualization and the encrypted one
		%final(1) = DC_Vizu+cryptDC;
%%End of the DC Encryption 

if(method==1)
		descramble = izigzag(final, sqrt(16), sqrt(16));
else 
		descramble = final;
end

end

