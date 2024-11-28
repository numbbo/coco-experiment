%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compiles the Coco functions in the Matlab/Octave wrapper  %
% to be able to run the exampleexperiment.m                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more off; % turn off page-wise output

fprintf("MEX     cocoCall.cpp\n");
mex("cocoCall.cpp");