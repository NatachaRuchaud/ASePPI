%---------------------------------------------------------
%% Compress and store the video frames -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the reference paper (Section 3): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs are:
% mov_yuv - video frames to compress and store
% Frame_start - which is the first frame
% ROI - the Region of interest (where the privacy protection is applied)
% imgCover - the cover that we use to keep the minimum information required
% by the video and to preserve the pleasantness of the video
% 
% outputs are:
% bitstream - the compressed data

function [bitstream] = compress_store(mov_yuv, Frame_start, ROI, imgCover, key, Quant, IntraPeriod, ext, bits, channel)
        Seq = mov_yuv(:,:,channel, Frame_start);
        Seq_Cover = imgCover(:,:,1);
        bitstream = '';
        %Save some needed values
        bitstream= bits;
        % Add '1111' to mark I frame 
        bitstream = [bitstream '1111'];
        %Encode I frame
        [Seq_i,bits_y]=encode_i_frame_demo_nat(Seq, Seq_Cover, Quant, ROI, channel, key);
        bitstream = [bitstream bits_y];
        X1 = Seq_i;
        %% Encode P frames  - subsequent frames
        for k = 2:IntraPeriod
            bitstream = [bitstream '0000'];
            % block_size = 16;        % Macroblock size for P frames
            [Seq_p,bits_y] = encode_p_frame_demo_nat(X1, mov_yuv(:,:,channel, k+Frame_start-1),Quant,ext,16, ROI, channel, key);
            bitstream = [bitstream bits_y];
            % the current reconstructed frame becomes the previous frame for
            % the next frame to be coded
            X1 = Seq_p;
        end 
end