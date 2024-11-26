%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compiles Matlab/Octave wrapper of the COCO functions      %
% to be able to run the exampleexperiment.m                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more off; % turn off page-wise output


fprintf('compiling cocoCall.c...');

mex("-s -Dchar16_t=uint16_t", "cocoCall.c");
mex("-s -Dchar16_t=uint16_t", "cocoProblemGetBestValue.c", "coco.c");

fprintf('Done\n');
fprintf('Preparation of all mex files finished.\n');