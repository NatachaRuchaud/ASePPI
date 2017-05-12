% *************************************************************************
%% H.264/AVC DECODER
% It works on YUV/AVI sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was first created by
% Abdullah Al Muhit
% contact - Almuhit [At] Gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M R Pickering, M R Frater and J F Arnold, �Video Coding using Elastic Motion Model and Larger Blocks,� IEEE Trans. Circ. And Syst. for Video Technology, vol. 20, no. 5, pp. 661-672, 2010. [Impact factor � 3.18] [PDF]
% A A Muhit, M R Pickering, M R Frater and J F Arnold, �Video Coding using Geometry Partitioning and an Elastic Motion Model,� accepted for publication in Journal of Visual Communication and Image Representation. [Impact factor � 1.33] [PDF]
% Modification is done by Ali Radmehr
% feel free to contact us: Radmehr.Telecom [AT] Gmail.Com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Taken back by Natacha Ruchaud to integrate her privacy protection method, ASePPI inside the H.264 framework
%% Please cite the following paper:  


%nameVideoOutput: the path of the coded video that you want to decode
%reverse: 1 if you want to recover the original video 0 otherwise
%(protected version)
%key: the same that you use in the encoding step mandatory if reverse = 1
%otherwise put a random value

%For example, decode_demo_nat('nat', 0, 0) or decode_demo_nat('grandma', 0, 0) to have the protected version
%decode_demo_nat('nat', 1, 123)  or decode_demo_nat('grandma', 1, 123) to reverse the process
%
%If face/people detection fails during the encoding the correct area may not protected


function decode_demo_nat(nameVideoOutput, reverse, key)
tic
system_dependent('DirChangeHandleWarn', 'Never');
addpath(genpath('.'));

global QP h w block_size

v = VideoWriter(['Videos_Resultats/', nameVideoOutput]); %, 'MPEG-4');  %Dufaux100Frames
open(v)
block_size = 16;        % Macroblock size for P frames

%% Input - load the encoded file
Frame_start = 1;
if(exist((['Res_Videos_Bitstream/', nameVideoOutput, '_bitstream_H264_enc', '.mat']))==2)
    bitstream_y = '';
    bitstream_u = '';
    bitstream_v = '';
    load(['Res_Videos_Bitstream/', nameVideoOutput, '_bitstream_H264_enc', '.mat']) % correct packets
    % Decode header
    [h,w,QP, IntraPeriod, Frame_start,LastFrame,m, ROI] = dec_header_demo(bitstream_y{Frame_start});
    IP = IntraPeriod;
    rgb_Ceq1 = [];
    rgb_Ceq2 = [];
    rgb_Ceq3 = [];
for iii = Frame_start:ceil(LastFrame/IntraPeriod)
    RoI = ROI((1+IntraPeriod*(iii-1)-1)/5+1, :);
      %---------------------------------------------------------
    rgb_Ceq = zeros(h, w, 3, IntraPeriod);
    
    r_k = zeros(h, w, IntraPeriod);
    if(LastFrame<iii*IntraPeriod)
        
        IntraPeriod =LastFrame-IntraPeriod*(iii-1)%(iii+IntraPeriod-Frame_end)
    end
    parfor channel=1:3
        idx = m;
        if(channel==1)
            bitstream = bitstream_y{iii};
            rgb_Ceq1 = [rgb_Ceq1; decompress(bitstream, idx, RoI, channel, reverse, key, IntraPeriod, h,w,QP)];
        elseif(channel==2)
            bitstream = bitstream_u{iii};
            rgb_Ceq2 = [rgb_Ceq2; decompress(bitstream, idx, RoI, channel, reverse, key, IntraPeriod, h,w,QP)];
        elseif(channel==3)
            bitstream = bitstream_v{iii};
            rgb_Ceq3 = [rgb_Ceq3; decompress(bitstream, idx, RoI, channel, reverse, key, IntraPeriod, h,w,QP)];
        end
        %rgb_Ceq2(1+IP*h*(iii-1):((iii-1)*IP+IntraPeriod)*h,:, channel) =decompress(bitstream, idx, RoI, channel, reverse, key, IntraPeriod, h,w,QP);
    end
end
    for k=Frame_start:LastFrame
        decompressedImg(:,:,1) = rgb_Ceq1(1+h*(k-1):h*(k),:,:);
        decompressedImg(:,:,2) = rgb_Ceq2(1+h*(k-1):h*(k),:,:);
        decompressedImg(:,:,3) = rgb_Ceq3(1+h*(k-1):h*(k),:,:);

        writeVideo(v, uint8(ycbcr2rgb(uint8(decompressedImg))));
        image(uint8(ycbcr2rgb(uint8(decompressedImg))))
        title(['Frame No. ' num2str(k)]);
        colormap(gray(256))
        truesize([2*h 2*w])
        drawnow
    end
toc
close(v)
end


end
