%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compiles the Coco functions in the Matlab/Octave wrapper  %
% to be able to run the exampleexperiment.m                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more off; % turn off page-wise output

fprintf("compiling cocoCall.c...");
mex("cocoCall.c", "coco.cpp");
fprintf('Done\n');
fprintf('Preparation of all mex files finished.\n');