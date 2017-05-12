%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Error Concealment methods
% Author: Ali Radmehr
% Contact me via : radmehr.telecom [AT] gmail.com
% Description: This function perform error concealment on erroneous frame
% which results in a better user experience in video quality.
% cur_frame = current frame
% prev_frame = previous frame
% mode = which method of error concealment should apply
% conceal_frame = the concealed frame using the mode give
% you can add other methods in other modes if needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conceal_frame = EC(cur_frame , prev_frame, ErrorMat,mode,block_size)

% controling...
if nargin~=5
    error('myApp:argChk', 'Wrong number of input arguments')
end

% default values
if isempty(mode)
    mode = 1;
end

switch mode
    
    case 1 % MB-copy method
        conceal_frame = MBCopy(cur_frame , prev_frame, ErrorMat,block_size);
end

function conceal_frame = MBCopy(cur_frame , prev_frame, ErrorMat,block_size)
% this function copy the mb in previous frame exactly in the same place
for i=1:length(ErrorMat(1,:))
    ii = ErrorMat(1,i);
    jj = ErrorMat(2,i);
    for k1 = (ii-1)*block_size+1:ii*block_size
        for k2 = (jj-1)*block_size+1:jj*block_size
            cur_frame(k1,k2)  = prev_frame(k1,k2);
        end
    end
end
conceal_frame = cur_frame;
