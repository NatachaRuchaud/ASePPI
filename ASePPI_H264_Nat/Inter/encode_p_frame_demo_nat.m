% Natacha Ruchaud added the SWRME (Search Window Restricted Motion Estimation) method
% Paper: Dai, F., Tong, L., Zhang, Y., & Li, J. (2011). Restricted H. 264/AVC video coding for privacy protected video scrambling. Journal of Visual Communication and Image Representation, 22(6), 479-490.

% inputs are:
% Seq1 - array containing the previous frame
% Seq2 - array containing the current frame to be coded
% Quant - the quantization parameter
% flag - flag to indicate whether to use extended cosine warping ME or not
% block_size - size of the macroBlock
% RoI - region of interest, the region to be protected
% Luma - If the channel is luminance 1 if chrominance 0
% key - The secret key to recover the original data with the decoder
%
% outputs are:
% Seq_r - the coded current frame
% bits_frame - the bits to save
% mse - the mean-squared-error of the frame
% R - the output bitrate for the frame

function [Seq_r,bits_frame] = encode_p_frame_demo_nat(Seq1,Seq2,Quant,flag, block_size, RoI, Luma, key)

% initialize some global variables:
% QP - the quantization parameter
% lambda_me - the lagrangian multiplier used in motion estimation
% lambda_mode - the lagrangian multiplier used in partition mode selection
% vlc_code_length - array containing code lengths for H.263 VLC tables
% vlc_code_length_last - array containing code lengths for H.263 VLC tables
% mvm - array containing previously calculated motion vectors
% S - array containing the previous frame and the current frame
% Sr - the coded current frame
% width - the width in pixels of the current frame
% height - the height in pixels of the current frame
% frame - the current frame to be coded
% extend - flag to indicate whether to use extended cosine warping ME
% is_ext - array used for results visualization

global QP lambda_me lambda_mode
global mvm S Sr width height frame extend is_ext

% H.264 global params
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros
load table.mat

% copy input variables to global variables (in MATLAB you cannot define global
% variables that have the same name as input variables)
Seq(:,:,1) = Seq1;
Seq(:,:,2) = Seq2;
S = Seq;
QP = Quant;
extend = flag;

% find the height, width and number of frames in the current input sequence
% (u is always 2 in this implementation)

[height,width,u] = size(S);

% calculate the lagrangian multipliers for the current value of QP

lambda_mode = 0.85*(2^((QP-12)/3));
lambda_me = lambda_mode^.5;


% initialize output variables

Sr = zeros(height,width,u);
Sr(:,:,1) = S(:,:,1);
mse = 0;
R = 0;
bits_frame = '';

% debugging and results visualization
% bits = zeros(height/16,width/16);
is_ext = zeros(height/4,width/4,4);


% figure(1)
% image(S(:,:,1)-S(:,:,2)+128)
% truesize([2*height 2*width])
% colormap(gray(256))
% drawnow

% in this implementation we are only coding P frames so "frame" is always 2

for frame = 2
    
    % initialize array for previously calulated motion vectors, infinity
    % indicates no motion vectors have been calculated for this pixel
    % horizontal vectors stored in mvm(:,:,1)
    % vertical vectors stored in mvm(:,:,2)
    
    mvm = inf.*ones(height+1,width+2,2);
    
    
    %     % debugging and results visualization
    %      figure(1)
    %      image(S(:,:,frame)-S(:,:,frame-1)+128)
    %      truesize([2*height 2*width])
    %      title('Encoding P Frame')
    %      drawnow
    
    % loop through the 16x16 macroblocks in the frame
    
    for i = 1:block_size:height
        for j = 1:block_size:width
            % bs - macroblock size
            % mr - motion vector search range
            
            %%%Force to be coded in 4*4 sub-macroblock partitions on the RoI to avoid drift error on the RoI border
            MB = false;
            if(Luma==1)
                if((i>RoI(2))&&(i<RoI(4))&&(j>RoI(1))&&(j<RoI(3)))
                    MB = true;
                end
            end
            bs = block_size;
            mr = 8;
            
            % MB header 0
            bits_frame = [bits_frame '0'];
            
            % determine the 16x16 partition mode
            % output variables:
            % rec_16 - the reconstructed 16x16 macroblock
            % R_16 - the output bitrate for this macroblock
            % ind_16 - partition index (1 = full, 2 = vert, 3 = horiz, 4 = quad)
            
            [rec_16,bits_16,ind_16,mvb] = mode_decision(bs,mr,i,j, RoI, Luma, MB, key);
            
            if ind_16 < 4   % if the partition mode is not quad
                
                % write the reconstructed 16x16 macroblock to the output
                
                Sr(i:i+block_size-1,j:j+block_size-1,frame) = rec_16;
                bits_frame = [bits_frame bits_16];
                
                % write the block motion vectors to the array containing the motion vectors
                % for the frame
                
                mvm(i:i+bs-1,j:j+bs-1,1) = mvb(:,:,1);
                mvm(i:i+bs-1,j:j+bs-1,2) = mvb(:,:,2);
                
            else    % else find the partition mode for each 8x8 sub-block
                
                % bs - macroblock size
                % mr - motion vector search range
                
                bs = block_size/2;
                mr = 4;
                % add header for a quad split
                bits_frame =[bits_frame '00100'];
                
                % determine the 8x8 partition mode for each sub-block,
                % write the reconstructed 8x8 sub-block to the output
                % array and add the rate for this sub-block to the rate
                % for the frame
                
                [rec_8,bits_8,ind_8,mvb] = mode_decision(bs,mr,i,j, RoI, Luma, MB, key);
                Sr(i:i+ bs -1,j:j+ bs -1,frame) = rec_8;;
                bits_frame = [bits_frame bits_8];
                
                mvm(i:i+bs-1,j:j+bs-1,1) = mvb(:,:,1);
                mvm(i:i+bs-1,j:j+bs-1,2) = mvb(:,:,2);
                
                
                [rec_8,bits_8,ind_8,mvb] = mode_decision(bs,mr,i,j+ bs, RoI, Luma, MB, key);
                Sr(i:i+ bs -1,j+ bs:j+ 2*bs -1,frame) = rec_8;
                bits_frame = [bits_frame bits_8];
                
                mvm(i:i+bs-1,j+bs:j+2*bs-1,1) = mvb(:,:,1);
                mvm(i:i+bs-1,j+bs:j+ 2*bs -1,2) = mvb(:,:,2);
                
                
                [rec_8,bits_8,ind_8,mvb] = mode_decision(bs,mr,i+ bs,j, RoI, Luma, MB, key);
                Sr(i+ bs:i+ 2*bs -1,j:j+ bs -1,frame) = rec_8;
                bits_frame = [bits_frame bits_8];
                mvm(i+ bs:i+ 2*bs -1,j:j+bs-1,1) = mvb(:,:,1);
                mvm(i+ bs:i+ 2*bs -1,j:j+bs-1,2) = mvb(:,:,2);
                
                [rec_8,bits_8,ind_8,mvb] = mode_decision(bs,mr,i+bs,j+bs, RoI, Luma, MB, key);
                Sr(i+ bs:i+ 2*bs -1,j+ bs:j+ 2*bs -1,frame) = rec_8;
                bits_frame = [bits_frame bits_8];
                
                mvm(i+ bs:i+ 2*bs -1,j+ bs:j+ 2*bs -1,1) = mvb(:,:,1);
                mvm(i+ bs:i+ 2*bs -1,j+ bs:j+ 2*bs -1,2) = mvb(:,:,2);
            end
        end
    end
    Seq_r = Sr(:,:,frame);
end



%--------------------------------------------------------------------------

function [rec,bits_b,ind,mvb] = mode_decision(bs,mr,i,j, RoI, Luma, MB, key)

global lambda_mode is_ext

% calculate the parameters for each partition mode
% output variables:
% mvb# - array containg motion vectors for the block
% rec# - array containing the reconstructed block
% mvr# - the bitrate required to transmit the motion vectors
% cfr# - the bitrate required to transmit the quantized coefficients
% ssd# - the sum-of-squared-difference for the block

[mvb0,rec0,bits0,ssd0] = motion_comp_full(bs,mr,i,j, RoI, Luma, key);

[mvb1,rec1,bits1,ssd1] = motion_comp_vert(bs,mr,i,j, RoI, Luma, key);

[mvb2,rec2,bits2,ssd2] = motion_comp_horz(bs,mr,i,j, RoI, Luma, key);

[mvb3,rec3,bits3,ssd3] = motion_comp_quad(bs,mr,i,j, RoI, Luma, key);

% mbtr# - the bitrate required to transmit the macroblock_type parameter

cfr0 = length(bits0);
cfr1 = length(bits1);
cfr2 = length(bits2);
cfr3 = length(bits3);

% find the partitioning mode which gives the minimum lagrangian cost
% function

[y,ind] = min([(ssd0+lambda_mode*cfr0) ...
    (ssd1+lambda_mode*cfr1) ...
    (ssd2+lambda_mode*cfr2) ...
    (ssd3+lambda_mode*cfr3)]);

% for the chosen partitioning mode write the partition variables to
% variables for the block
if(MB)
    ind = 4;
end
switch ind
    case 1
        rec = rec0;
        bits_b = bits0;
        ssd = ssd0;
        mvb = mvb0;
    case 2
        rec = rec1;
        bits_b = bits1;
        ssd = ssd1;
        mvb = mvb1;
    case 3
        rec = rec2;
        bits_b = bits2;
        ssd = ssd2;
        mvb = mvb2;
    case 4
        rec = rec3;
        bits_b = bits3;
        ssd = ssd3;
        mvb = mvb3;
end



%--------------------------------------------------------------------------

function [mvb,rec,bits_b,ssd] = motion_comp_full(bs,mr,i,j, RoI, Luma, key)

global mode

% mode - variable indicating which partitioning mode to use for the block

mode = 1;

bits_b = '1';

% calculate the motion parameters for one partition consisting of the whole
% block
% mvx - horizontal motion vector
% mvy - vertical motion vector

[mvx,mvy,rec,bits,ssd] = motion_comp_block(bs,bs,mr,i,j,'m', RoI, Luma, key);
bits_b = [bits_b bits];

mvb(:,:,1) = mvx.*ones(bs,bs);
mvb(:,:,2) = mvy.*ones(bs,bs);


%--------------------------------------------------------------------------

function [mvb,rec,bits_b,ssd] = motion_comp_vert(bs,mr,i,j, RoI, Luma, key)

global mode

% mode - variable indicating which partitioning mode to use for the block

mode = 2;
bits_b = '010';

% calculate the motion parameters for two vertical partitions
% mvx# - horizontal motion vector for partition #
% mvy# - vertical motion vector for partition #

[mvx0,mvy0,rec0,bits0,ssd0] = motion_comp_block(bs/2,bs,mr,i,j,'a', RoI, Luma, key);
[mvx1,mvy1,rec1,bits1,ssd1] = motion_comp_block(bs/2,bs,mr,i,j+bs/2,'c', RoI, Luma, key);

% calculate the combined parameters for the whole block

rec = [rec0 rec1];
bits_b = [bits_b bits0 bits1];
ssd = ssd0+ssd1;

mvb(:,:,1) = [mvx0.*ones(bs,bs/2) mvx1.*ones(bs,bs/2)];
mvb(:,:,2) = [mvy0.*ones(bs,bs/2) mvy1.*ones(bs,bs/2)];


%--------------------------------------------------------------------------

function [mvb,rec,bits_b,ssd] = motion_comp_horz(bs,mr,i,j, RoI, Luma, key)

global mode

% mode - variable indicating which partitioning mode to use for the block

mode = 3;
bits_b = '011';

% calculate the motion parameters for two horizontal partitions
% mvx# - horizontal motion vector for partition #
% mvy# - vertical motion vector for partition #

[mvx0,mvy0,rec0,bits0,ssd0] = motion_comp_block(bs,bs/2,mr,i,j,'b', RoI, Luma, key);
[mvx1,mvy1,rec1,bits1,ssd1] = motion_comp_block(bs,bs/2,mr,i+bs/2,j,'a', RoI, Luma, key);

% calculate the combined parameters for the whole block

rec = [rec0;rec1];
bits_b = [bits_b bits0 bits1];
ssd = ssd0+ssd1;

mvb(:,:,1) = [mvx0.*ones(bs/2,bs); mvx1.*ones(bs/2,bs)];
mvb(:,:,2) = [mvy0.*ones(bs/2,bs); mvy1.*ones(bs/2,bs)];


%--------------------------------------------------------------------------

function [mvb,rec,bits_b,ssd] = motion_comp_quad(bs,mr,i,j, RoI, Luma, key)

global mode

% mode - variable indicating which partitioning mode to use for the block

mode = 4;
bits_b = '00100';

% calculate the motion parameters for four partitions
% mvx# - horizontal motion vector for partition #
% mvy# - vertical motion vector for partition #

[mvx0,mvy0,rec0,bits0,ssd0] = motion_comp_block(bs/2,bs/2,mr,i,j,'m', RoI, Luma, key);

[mvx1,mvy1,rec1,bits1,ssd1] = motion_comp_block(bs/2,bs/2,mr,i,j+bs/2,'m', RoI, Luma, key);

[mvx2,mvy2,rec2,bits2,ssd2] = motion_comp_block(bs/2,bs/2,mr,i+bs/2,j,'m', RoI, Luma, key);

[mvx3,mvy3,rec3,bits3,ssd3] = motion_comp_block(bs/2,bs/2,mr,i+bs/2,j+bs/2,'m', RoI, Luma, key);

% calculate the combined parameters for the whole block

rec = [rec0 rec1;rec2 rec3];

bits_b = [bits_b bits0 bits1 bits2 bits3];
ssd = ssd0+ssd1+ssd2+ssd3;

mvb(:,:,1) = [mvx0.*ones(bs/2,bs/2) mvx1.*ones(bs/2,bs/2);mvx2.*ones(bs/2,bs/2) mvx3.*ones(bs/2,bs/2)];
mvb(:,:,2) = [mvy0.*ones(bs/2,bs/2) mvy1.*ones(bs/2,bs/2);mvy2.*ones(bs/2,bs/2) mvy3.*ones(bs/2,bs/2)];



%--------------------------------------------------------------------------

function [mvx,mvy,rec,bits_b,ssd] = motion_comp_block(bsx,bsy,mr,i,j,mvp_mode, RoI, Luma, key);

global lambda_mode mvm S Sr width height frame extend is_ext mode

bits_b = '0';

% extract the search window from the previous reconstructed frame
swr = mirror_pad(extract_object(Sr(:,:,frame-1),[j-mr i-mr bsx+2*mr bsy+2*mr],[1 1 width height], 0));

% extract the current block from the previous reconstructed frame
curr = extract_object(S(:,:,frame),[j i bsx bsy],[1 1 width height], 0);

% calculate the motion vector predictions mvxp andmvyp
switch mvp_mode
    case 'm'
        [mvxp,mvyp] = find_median_mv_pred(i,j,bsx);
    case 'a'
        [mvxp,mvyp] = find_mv_pred_a(i,j,bsx);
    case 'b'
        [mvxp,mvyp] = find_mv_pred_b(i,j,bsx);
    case 'c'
        [mvxp,mvyp] = find_mv_pred_c(i,j,bsx);
end

% calculate the motion vectors and the motion compensated prediction for
% the block
% mvx - horizontal motion vector
% mvy - vertical motion vector
% pred - motion compensated prediction block

if(Luma)
    [mvx,mvy,pred] = motion_comp_Luma(swr,curr,bsx,bsy,mr,mvxp,mvyp, RoI, i, j);
else
    [mvx,mvy,pred] = motion_comp(swr,curr,bsx,bsy,mr,mvxp,mvyp, RoI, i, j);
end

% updated the frame motion vector arrays so that motion vectors from
% previously coded partitions of a block can be used as predictions

% mvm(i+1:i+bsy,j+1:j+bsx,1) = mvx.*ones(bsy,bsx);
% mvm(i+1:i+bsy,j+1:j+bsx,2) = mvy.*ones(bsy,bsx);

% calculate the bitrate required to transmit the motion vectors
% mvr = code_mv(mvx-mvxp,mvy-mvyp);
bits_b = [bits_b enc_golomb(4*(mvx-mvxp),1)];
bits_b = [bits_b enc_golomb(4*(mvy-mvyp),1)];

% calculate the motion-compensated difference block mcp

mcp = curr-pred;


% calculate the bitrate required to transmit the motion-compensated
% prediction block cfr and the reconstructed motion-compensated
% prediction block mcp_r

[mcp_r,bits] = code_block(mcp);

%%If luminance channel and inside the RoI use ASePPI (i.e. the privacy protection tool)
if(Luma==1)
    if((i>RoI(2))&&(i<RoI(4))&&(j>RoI(1))&&(j<RoI(3)))
        [mcp_r2,bits] = code_block_protect(mcp, Luma, RoI, i, j, key);
    end
end
bits_b = [bits_b bits];

% calculate the reconstucted block rec
rec = mcp_r+pred;


% calculate the sum-of-squared-difference for the reconstructed block
ssd = round(sum(sum((curr-rec).^2)));
sad = sum(sum(abs(curr-rec)));

% mvr_n - the bitrate required for the motion vectors for normal ME
% cfr_n - the bitrate required for the quantized coefficients for normal ME

% mvr_n = mvr;
% cfr_n = cfr;



%--------------------------------------------------------------------------

function [mvx,mvy,pred] = motion_comp(swr,curr,bsx,bsy,mr,mvxp,mvyp, RoI, ii, jj)

global lambda_me sw_sup swr_sup


min_J = inf;

[n,m] = size(swr);
[nc,mc] = size(curr);


% Create super-resolution version of the original and coded search windows
% by interpolating each by a factor of 4. The function mirror_pad
% symmetrically extends data from inside the frame for those search windows
% which extend outside the frame boundary.

% sw_sup = interp2([1:m],[1:n]',mirror_pad(sw),[1:.25:m],[1:.25:n]','bicubic');
swr_sup = interp2([1:m],[1:n]',mirror_pad(swr),[1:.25:m],[1:.25:n]','bicubic');


[n,m] = size(swr_sup);

% loop through each integer pixel location at the original resolution

for i = 5:4:n-4*bsy
    for j = 5:4:m-4*bsx
        
        % calculate the number of bits required for motion vectors at this
        % candidate location
        
        mvr = code_mv((j-1)/4-mr-mvxp,(i-1)/4-mr-mvyp);
        
        
        % calculate the sum-of-absolute-difference between the current
        % block and the block at this candidate location
        
        sad = sum(sum(abs(swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1)-curr)));
        
        % calulate the lagrangian cost function for this location
        
        J = sad+lambda_me*mvr;
        
        % if the lagrangian cost function is minimised at this location
        % assign the current motion vectors to the integer pixel motion
        % vectors for this block
        
        if J < min_J
            min_J = J;
            imvx = j;
            imvy = i;
        end
    end
end

% Full pixel accuracy only
j=imvx;
i=imvy;
mvx = (j-1)/4-mr;
mvy = (i-1)/4-mr;
pred = swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1);
% full pixel

min_J = inf;

% loop through each quarter pixel location around the integer pixel
% location from the previous step

for i = imvy-3:imvy+3
    for j = imvx-3:imvx+3
        
        
        % calculate the number of bits required for motion vectors at this
        % candidate location
        
        mvr = code_mv((j-1)/4-mr-mvxp,(i-1)/4-mr-mvyp);
        
        
        % calculate the sum-of-absolute-difference between the current
        % block and the block at this candidate location
        
        sad = sum(sum(abs(swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1)-curr)));
        
        
        % calulate the lagrangian cost function for this location
        
        J = sad+lambda_me*mvr;
        
        
        % if the lagrangian cost function is minimised at this location
        % assign the current motion vectors to the quarter pixel motion
        % vectors for this block and assign the interpolated block at this
        % location to the motion compensated prediction for this block
        
        if J < min_J
            min_J = J;
            mvx = (j-1)/4-mr;
            mvy = (i-1)/4-mr;
            pred = swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1);
        end
        
    end
end





%--------------------------------------------------------------------------

function [mvx,mvy,pred] = motion_comp_Luma(swr,curr,bsx,bsy,mr,mvxp,mvyp, RoI, ii, jj)

global lambda_me sw_sup swr_sup

min_J = inf;

[n,m] = size(swr);
[nc,mc] = size(curr);

% Create super-resolution version of the original and coded search windows
% by interpolating each by a factor of 4. The function mirror_pad
% symmetrically extends data from inside the frame for those search windows
% which extend outside the frame boundary.
% sw_sup = interp2([1:m],[1:n]',mirror_pad(sw),[1:.25:m],[1:.25:n]','bicubic');
swr_sup = interp2([1:m],[1:n]',mirror_pad(swr),[1:.25:m],[1:.25:n]','bicubic');
[n,m] = size(swr_sup);


%Integrate the SWRME method: the inter prediction from the privacy region to the non-privacy region has forbidden

lc = jj;
rc = jj+bsx;
tc = ii;
bc = ii+bsy;
RoI = RoI+1;

%%%Outside RoI
lcOut = ((lc<RoI(1))||(lc>=RoI(3)));
rcOut = ((rc<=RoI(1))||(rc>RoI(3)));
tcOut = ((tc<RoI(2))||(tc>=RoI(4)));
bcOut = ((bc<=RoI(2))||(bc>RoI(4)));

% loop through each integer pixel location at the original resolution
for i = 5:4:n-4*bsy
    for j = 5:4:m-4*bsx
        %Current block is Outside RoI
        if((lcOut&&rcOut)||(tcOut&&bcOut))
            lp = jj-mr+(j-1)/4-1;
            rp = lp+bsx;
            tp = ii-mr+(i-1)/4-1;
            bp = tp+bsy;
            lpOut = ((lp<RoI(1)-2)||(lp>=RoI(3)));
            rpOut = ((rp<=RoI(1)-2)||(rp>RoI(3)));
            tpOut = ((tp<RoI(2)-2)||(tp>=RoI(4)));
            bpOut = ((bp<=RoI(2)-2)||(bp>RoI(4)));
            
            %%%%%Faire only IF previous block is Totally outside RoI !!!
            if((lpOut&&rpOut)||(tpOut&&bpOut))
                % calculate the number of bits required for motion vectors at this
                % candidate location
                mvr = code_mv((j-1)/4-mr-mvxp,(i-1)/4-mr-mvyp);
                
                % calculate the sum-of-absolute-difference between the current
                % block and the block at this candidate location
                sad = sum(sum(abs(swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1)-curr)));
                
                % calulate the lagrangian cost function for this location
                J = sad+lambda_me*mvr;
                
                % if the lagrangian cost function is minimised at this location
                % assign the current motion vectors to the integer pixel motion
                % vectors for this block
                if J < min_J
                    min_J = J;
                    imvx = j;
                    imvy = i;
                end
                
            end
            %%IF current block is inside RoI
        else
            % calculate the number of bits required for motion vectors at this
            % candidate location
            
            mvr = code_mv((j-1)/4-mr-mvxp,(i-1)/4-mr-mvyp);
            
            % calculate the sum-of-absolute-difference between the current
            % block and the block at this candidate location
            
            sad = sum(sum(abs(swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1)-curr)));
            
            % calulate the lagrangian cost function for this location
            
            J = sad+lambda_me*mvr;
            
            % if the lagrangian cost function is minimised at this location
            % assign the current motion vectors to the integer pixel motion
            % vectors for this block
            
            if J < min_J
                min_J = J;
                imvx = j;
                imvy = i;
            end
        end
    end
end
% Full pixel accuracy only
j=imvx;
i=imvy;

mvx = (j-1)/4-mr;
mvy = (i-1)/4-mr;

pred = swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1);
% full pixel

min_J = inf;

% loop through each quarter pixel location around the integer pixel
% location from the previous step
for i = imvy-3:imvy+3
    for j = imvx-3:imvx+3
        %Current block is Outside RoI
        if((lcOut&&rcOut)||(tcOut&&bcOut))
            lp = jj-mr+(j-1)/4-1;
            rp = lp+bsx;
            tp = ii-mr+(i-1)/4-1;
            bp = tp+bsy;%
            lpOut = ((lp<RoI(1)-2)||(lp>=RoI(3)));
            rpOut = ((rp<=RoI(1)-2)||(rp>RoI(3)));
            tpOut = ((tp<RoI(2)-2)||(tp>=RoI(4)));
            bpOut = ((bp<=RoI(2)-2)||(bp>RoI(4)));
            
            %%%%%Faire only IF previous block is outside RoI !!!
            if((lpOut&&rpOut)||(tpOut&&bpOut))
                
                % calculate the number of bits required for motion vectors at this
                % candidate location
                mvr = code_mv((j-1)/4-mr-mvxp,(i-1)/4-mr-mvyp);
                % calculate the sum-of-absolute-difference between the current
                % block and the block at this candidate location
                sad = sum(sum(abs(swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1)-curr)));
                % calulate the lagrangian cost function for this location
                J = sad+lambda_me*mvr;
                % if the lagrangian cost function is minimised at this location
                % assign the current motion vectors to the quarter pixel motion
                % vectors for this block and assign the interpolated block at this
                % location to the motion compensated prediction for this block
                if J < min_J
                    min_J = J;
                    mvx = (j-1)/4-mr;
                    mvy = (i-1)/4-mr;
                    pred = swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1);
                end
            end
            %%IF current block is inside RoI
        else
            % calculate the number of bits required for motion vectors at this
            % candidate location
            mvr = code_mv((j-1)/4-mr-mvxp,(i-1)/4-mr-mvyp);
            % calculate the sum-of-absolute-difference between the current
            % block and the block at this candidate location
            
            sad = sum(sum(abs(swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1)-curr)));
            
            % calulate the lagrangian cost function for this location
            
            J = sad+lambda_me*mvr;
            
            % if the lagrangian cost function is minimised at this location
            % assign the current motion vectors to the quarter pixel motion
            % vectors for this block and assign the interpolated block at this
            % location to the motion compensated prediction for this block
            
            if J < min_J
                min_J = J;
                mvx = (j-1)/4-mr;
                mvy = (i-1)/4-mr;
                pred = swr_sup(i:4:i+4*bsy-1,j:4:j+4*bsx-1);
            end
        end
        
    end
end


%--------------------------------------------------------------------------
%% Transform, Quantization, Entropy by Muhit
% transform = Integer transform
% Quantization = h.264
% VLC = CAVLC (H.264)

function [err_r,bits_b] = code_block(err)

global QP

[n,m] = size(err);
bits_b = '';
cfr = 0;

for i = 1:4:n
    for j = 1:4:m
        
        c(i:i+3,j:j+3) = integer_transform(err(i:i+3,j:j+3));
        cq(i:i+3,j:j+3) = quantization(c(i:i+3,j:j+3),QP);
        [bits] = enc_cavlc(cq(i:i+3,j:j+3), 0, 0);
        bits_b = [bits_b bits];
        
        Wi = inv_quantization(cq(i:i+3,j:j+3),QP);
        
        Y = inv_integer_transform(Wi);
        
        err_r(i:i+3,j:j+3) = round(Y/64);
        
    end
end

function [err_r,bits_b] = code_block_protect(err, luma, RoI, ii, jj, key)

global QP

[n,m] = size(err);
bits_b = '';
cfr = 0;
seed = ii+key;
for i = 1:4:n
    for j = 1:4:m
        c(i:i+3,j:j+3) = integer_transform(err(i:i+3,j:j+3));
        Size = size(c(i:i+3,j:j+3));
        cq(i:i+3,j:j+3) = quantization(c(i:i+3,j:j+3),QP);
        
        block = cq(i:i+3,j:j+3);
        %%Scrambling step
        [scramble] = scrambledXoR_Inter_Forward_demo_nat(block, seed);
        cq(i:i+3,j:j+3) = scramble;
        
        [bits] = enc_cavlc(cq(i:i+3,j:j+3), 0, 0);
        bits_b = [bits_b bits];
        
        Wi = inv_quantization(cq(i:i+3,j:j+3),QP);
        Y = inv_integer_transform(Wi);
        err_r(i:i+3,j:j+3) = round(Y/64);
    end
end

%--------------------------------------------------------------------------

function mvr = code_mv(mvx,mvy)

if mvx < 0
    code_num = 8*abs(mvx);
elseif mvx > 0
    code_num = 8*abs(mvx)-1;
else
    code_num = 0;
end

M = floor(log2(code_num+1));

mvr = 2*M+1;

if mvy < 0
    code_num = 8*abs(mvy);
elseif mvy > 0
    code_num = 8*abs(mvy)-1;
else
    code_num = 0;
end

M = floor(log2(code_num+1));

mvr = mvr+2*M+1;


%--------------------------------------------------------------------------

function mvr = code_ext_mv(mvx,mvy,p)

if mvx < 0
    code_num = 2*p*abs(mvx);
elseif mvx > 0
    code_num = 2*p*abs(mvx)-1;
else
    code_num = 0;
end

M = floor(log2(code_num+1));

mvr = 2*M+1;

if mvy < 0
    code_num = 2*p*abs(mvy);
elseif mvy > 0
    code_num = 2*p*abs(mvy)-1;
else
    code_num = 0;
end

M = floor(log2(code_num+1));

mvr = mvr+2*M+1;

%--------------------------------------------------------------------------

function mbtr = code_mb_type(mb_type);

M = floor(log2(code_num+1));

mbtr = 2*M+1;



%--------------------------------------------------------------------------

function [mvxp,mvyp] = find_median_mv_pred(i,j,bs);

global mvm

for k = 1:2
    mvpa = mvm(i+1,j,k);
    mvpb = mvm(i,j+1,k);
    mvpc = mvm(i,j+bs+1,k);
    mvpd = mvm(i,j,k);
    
    % bin2dec([num2str(isinf(mvpa)) num2str(isinf(mvpb)) num2str(isinf(mvpc)) num2str(isinf(mvpd))])
    switch bin2dec([num2str(isinf(mvpa)) num2str(isinf(mvpb)) num2str(isinf(mvpc)) num2str(isinf(mvpd))])
        case 0
            mvp(k) = median([mvpa mvpb mvpc]);
        case 1
            mvp(k) = median([mvpa mvpb mvpc]);
        case 2
            mvp(k) = median([mvpa mvpb mvpd]);
        case 3
            mvp(k) = median([mvpa mvpb]);
        case 4
            mvp(k) = median([mvpa mvpc mvpd]);
        case 5
            mvp(k) = median([mvpa mvpc]);
        case 6
            mvp(k) = median([mvpa mvpd]);
        case 7
            mvp(k) = mvpa;
        case 8
            mvp(k) = median([mvpb mvpc mvpd]);
        case 9
            mvp(k) = median([mvpb mvpc]);
        case 10
            mvp(k) = median([mvpb mvpd]);
        case 11
            mvp(k) = mvpb;
        case 12
            mvp(k) = median([mvpc mvpd]);
        case 13
            mvp(k) = mvpc;
        case 14
            mvp(k) = mvpd;
        case 15
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);
%--------------------------------------------------------------------------

function [mvxp,mvyp] = find_mv_pred_a(i,j,bs);

global mvm

for k = 1:2
    
    mvpa = mvm(i+1,j,k);
    
    switch isinf(mvpa)
        case 0
            mvp(k) = mvpa;
        case 1
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);



%--------------------------------------------------------------------------

function [mvxp,mvyp] = find_mv_pred_b(i,j,bs);

global mvm

for k = 1:2
    
    mvpb = mvm(i,j+1,k);
    
    switch isinf(mvpb)
        case 0
            mvp(k) = mvpb;
        case 1
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);



%--------------------------------------------------------------------------

function [mvxp,mvyp] = find_mv_pred_c(i,j,bs);

global mvm

for k = 1:2
    
    mvpc = mvm(i,j+bs+1,k);
    
    switch isinf(mvpc)
        case 0
            mvp(k) = mvpc;
        case 1
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);



%--------------------------------------------------------------------------



