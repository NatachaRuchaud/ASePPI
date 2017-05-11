# ASePPI
Adaptive Scrambling enabling Privacy Protection and Intelligibility 

The code integrates a privacy protection method, ASePPI, into the H.264/AVC codec and runs on Matlab. 


To encode a video (.avi or .yuv) use the function encode_demo_nat(nameVideo, face, nameVideoOutput, key).
In the folder "videos/", you will find several video examples to encode.

For example, do:
encode_demo_nat('videos/grandma_qcif.yuv', 1, 'grandma', 123)



To decode a video (that you have encoded) use the function decode_demo_nat(nameVideoOutput, reverse, key)
In the folder "Res_Videos_Bitstream/", you will find the videos that you have encoded.

For example, do:
decode_demo_nat('grandma', 0, 0) to have the protected version
decode_demo_nat('grandma', 1, 123) to reverse the process and almost obtain the original video
