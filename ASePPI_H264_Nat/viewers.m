%---------------------------------------------------------
%% Watching the privacy protected video in real time -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by
% Natacha Ruchaud
% contact - ruchaud@eurecom.fr
% website - https://eurecom.fr/~ruchaud
% Please cite the reference paper (Section 3.2.4): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inputs are:
% mov - video frames to protect
% 
% outputs are:
% mov_yuv - original video frames to protect
% RoI - the Region of interest (where the privacy protection is applied)
% imgCover - the cover that we use to keep the minimum information required


function [mov_yuv, imgCover, RoI] = viewers(mov, face, detector, peopleDetector, IntraPeriod, Frame_start, LastFrame)
RoI = [];
for i =Frame_start:LastFrame
        mov_yuv(:,:,:,i) = double(mov(i+0).cdata);
        
        %%%To reduce the image size by 4
        %[hh, ww, cc]=size(mov(i+0).cdata);
        %hh = round(hh/64)*16;
        %ww = round(ww/64)*16;
        %mov_yuv(:,:,:,i) = imresize(double(mov(i+0).cdata), [hh, ww]);
        %%%To reduce the image size by 4
        
        %Apply the detector on the orignal size to have better resutls
        if(face)
            BBOXES = step(detector,ycbcr2rgb(mov(i+0).cdata));
        else
            BBOXES = step(peopleDetector,ycbcr2rgb(mov(i+0).cdata));
        end
        %%%If the image size is changed
        %BBOXES = BBOXES/4;

    if(mod(i, IntraPeriod)==0)||(i==LastFrame)
        ROI = [min(BBOXES(:, 1)), min(BBOXES(:, 2)), max(BBOXES(:, 3))+min(BBOXES(:,1)), max(BBOXES(:, 4))+min(BBOXES(:,2))];
        ROI =  round(ROI/4)*4;
        if(ROI(1)==0)
            ROI(1)=4;
        end
        if(ROI(2)==0)
            ROI(2)=4;
        end
        RoI = [RoI; ROI];
        % Compute the size of the pixelization for the appearance of the final
        % image
        height_RoI = ROI(4)-ROI(2);
        width_RoI =  ROI(3)-ROI(1);
        kn = ceil(sqrt((width_RoI*height_RoI)/99)/4)*4;
        km = kn;
        for ii =i-IntraPeriod+1:i
            %%What the viewers see directly
            imgCover(:,:,:, ii) =pix(mov_yuv(:,:,:,ii), kn, km);
            imgCoverRoI = imgCover(ROI(2)+1:ROI(4), ROI(1)+1:ROI(3), :, ii);
            imgCover(:,:,:, ii) = mov_yuv(:,:,:,ii);
            imgCover(ROI(2)+1:ROI(4), ROI(1)+1:ROI(3), 1, ii) =imgCoverRoI(:,:,1);
            imshow(uint8(ycbcr2rgb(uint8(imgCover(:,:,:, ii)))))
            title(['Frame No. ' num2str(ii)]);
            colormap(gray(256))
            drawnow
        end
    end
end
end