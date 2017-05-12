% *************************************************************************
%% H.264/AVC ENCODER
% It works on YUV/AVI sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was first created by
% Abdullah Al Muhit
% contact - Almuhit [At] Gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M R Pickering, M R Frater and J F Arnold, Video Coding using Elastic Motion Model and Larger Blocks, IEEE Trans. Circ. And Syst. for Video Technology, vol. 20, no. 5, pp. 661-672, 2010. [Impact factor � 3.18] [PDF]
% A A Muhit, M R Pickering, M R Frater and J F Arnold, Video Coding using Geometry Partitioning and an Elastic Motion Model, accepted for publication in Journal of Visual Communication and Image Representation. [Impact factor � 1.33] [PDF]
% Modification is done by Ali Radmehr
% feel free to contact us: Radmehr.Telecom [AT] Gmail.Com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Taken back/Improved by Natacha Ruchaud to integrate her privacy protection method, ASePPI inside the H.264 framework
% It works on polychrome image sequence
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the following paper: 


%nameVideo: the path of your video that you want to test (.avi)
%face: if the video contain a face 1, a body 0
%nameVideoOutput: the name of the encoded video Output
%key: the secret key that you need to recover the original video, feel free
%to choose the number that you want

%For example, run encode_demo_nat('videos/Wave_arms_Original_Video.avi', 0, 'nat', 123)
%or encode_demo_nat('videos/grandma_qcif.yuv', 1, 'grandma', 123)

%You can download rows videos (.yuv) in QCIF at
%http://trace.eas.asu.edu/yuv/,

%If the detection of face or people fails our process will not protect the
%correct area

function encode_demo_nat(nameVideo, face, nameVideoOutput, key)
tic
addpath(genpath('.'));
%Face detector
detector = vision.CascadeObjectDetector('ClassificationModel', 'FrontalFaceLBP', 'MergeThreshold', 0);
%People Detector
peopleDetector = vision.PeopleDetector('ClassificationThreshold', 0);

addpath('readVideos/');

Quant =12; % Quantization parameter
IntraPeriod =5; % Intra Period
ext = 1;                % switch for extended M
Frame_start=1;

if(strcmp(nameVideo(end-2:end),'avi'))
    [mov] = loadFileAvi(nameVideo); % loading avi file
    LastFrame = length(mov);
elseif(strcmp(nameVideo(end-2:end),'yuv'))%QCIF size
    width = 176;%352;
    height = 144;%288;
    LastFrame=100;
    [mov] = loadFileYuv(nameVideo,width,height,Frame_start:LastFrame); % loading yuv file
end

%%%Detection of the privacy area (Body or Face) every IP
%%%Create the cover images and show a real time protected video version
[mov_yuv, imgCover, RoI] = viewers(mov, face, detector, peopleDetector, IntraPeriod, Frame_start, LastFrame);

[h,w,u] = size(mov_yuv(:,:,1,1));
%Save some needed values
[bits] = header_demo(h,w,Quant, IntraPeriod, Frame_start, LastFrame, RoI);

%%%Compress and store the data
for iii = Frame_start:ceil(LastFrame/IntraPeriod)
    if(LastFrame<iii*IntraPeriod)
        IntraPeriodModify =LastFrame-IntraPeriod*(iii-1)
    else
        IntraPeriodModify = IntraPeriod;
    end
    %---------------------------------------------------------  
    bitstream_yy = [];
    bitstream_uu = [];
    bitstream_vv = [];
    parfor channel=1:3
        [bitstream] = compress_store(mov_yuv, iii, RoI(iii, :), imgCover(:,:,:, 1+IntraPeriod*(iii-1)), key, Quant, IntraPeriodModify, ext, bits, channel);
        if(channel==1)
            bitstream_yy = [bitstream_yy; bitstream];
        elseif(channel==2)
            bitstream_uu = [bitstream_uu; bitstream];
        elseif(channel==3)
            bitstream_vv = [bitstream_vv; bitstream];
        end
    end
    bitstream_y{iii} = bitstream_yy;
    bitstream_u{iii} = bitstream_uu;
    bitstream_v{iii} = bitstream_vv;
end
%Save/Store the Output
save(['Res_Videos_Bitstream/', nameVideoOutput, '_bitstream_H264','_enc'], 'bitstream_y', 'bitstream_u', 'bitstream_v')
toc

