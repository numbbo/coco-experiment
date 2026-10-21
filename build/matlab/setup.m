%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compiles the Coco functions in the Matlab/Octave wrapper  %
% to be able to run the exampleexperiment.m                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more off; % turn off page-wise output

% The COCO library is C++ and is compiled separately because the char16_t
% define needed for cocoCall.c on some systems must not be applied to C++.
% coco.cpp needs C++17, which is not the default of every compiler that mex uses.
fprintf('compiling coco.cpp...');
if exist('OCTAVE_VERSION', 'builtin')
  setenv('CXXFLAGS', [strtrim(mkoctfile('-p', 'CXXFLAGS')) ' -std=c++17']);
  mex('-c', 'coco.cpp');
elseif ispc
  mex('COMPFLAGS=$COMPFLAGS /std:c++17', '-c', 'coco.cpp');
else
  mex('CXXFLAGS=$CXXFLAGS -std=c++17', '-c', 'coco.cpp');
end
fprintf('Done\n');

fprintf('compiling cocoCall.c...');
coco_object = 'coco.o';
if ~exist(coco_object, 'file')
  coco_object = 'coco.obj'; % MSVC
end
mex('-Dchar16_t=uint16_t', 'cocoCall.c', coco_object);
fprintf('Done\n');
fprintf('Preparation of all mex files finished.\n');