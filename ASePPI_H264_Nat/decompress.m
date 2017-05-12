%---------------------------------------------------------
%% Decompress and show the video frames -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% The Paper ASePPI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs are:
% bitstream - the compressed video frames
% idx - index to read the bitstream
% RoI - the Region of interest (where the privacy protection is applied)
% 
% outputs are:
% bitstream - the compressed data


function [decompressed_data] = decompress(bitstream, idx, RoI, channel, reverse, key, IntraPeriod, h, w, QP)
        if (strcmp(bitstream(idx:idx+3),'1111'))
            disp('Decoding I Frame y')
            idx = idx + 4;
            [Ceq(:,:,channel, 1),idx]=decode_i_frame_demo_nat(idx,bitstream, RoI, channel, reverse, key, h,w,QP);
        end
        X1(:,:,channel) = Ceq(:,:,channel, 1);
        
        R =[];
        R = [R; X1(:,:,channel)];
        
        for k = 2:IntraPeriod
            if (strcmp(bitstream(idx:idx+3),'0000'))
                disp('Decoding P Frame y')
                idx = idx + 4;
                %-----------PERFORMING EC ON THE ERRONEOUS FRAME---------------
                B=0; %1.4; % Average Burst Length
                PLR=0; %5/100; % Packet Loss Rate
                [r,idx,ErrorMat]= decode_p_frame_demo_nat(idx,bitstream,X1(:,:,channel),B,PLR, RoI, channel, k, reverse, key); % decoding using ec
                if isempty(ErrorMat) == 0 % if there is any error, EC would be applied.
                    r = EC(r, X1(:,:,channel),ErrorMat,1,16); %mode 1 is frame-copy method.
                end
                R =[R; r];
                X1(:,:,channel) = r;
            end
        end
        decompressed_data = R;
end