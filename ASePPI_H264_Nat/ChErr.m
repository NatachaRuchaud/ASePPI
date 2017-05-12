% The code below simulating a channel to make erroneous packets for
% H.264/AVC video codec

function states = ChErr( B, PLR)
%bits_err_less = bitstream;
%load ./bitstream_enc18_11.mat % load your bitsream enc output! <!-- it should be loaded -->

% ------Dr. Ghanbari Simulation Code for Genrating Errors-------
% 99 states which have markov model, according to his marvelous book
%B=1.4; % Average Burst Length
%PLR= 5/100; % Packet Loss Rate
Alpha= 1-PLR; %1/B;
Beta= PLR/(B*(1-PLR));
PreviousPacketLost=0; % 0=loss
%states=zeros(1,99);
%for i=1: 4 :length(bitstream)-4 % for 99 MBs in a frame. for qcif nd 16x16 MB
    switch PreviousPacketLost
        case 0
            if rand< (1-Alpha)
                states=0;
                %bitstream(i)=0;
            else
                states=1;
            end
        case 1
            if rand< Beta
                states=0;
                %bitstream(i)=0;
            else
                states=1;
            end
            
    end
    PreviousPacketLost=states;
%end

%save(['bitstream','_enc','_err'],'bitstream')



