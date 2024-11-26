% Returns the best object value that can be achieved.
%
% Parameters:
%   problem  The given problem.
function value = cocoProblemGetBestValue(problem)
value = cocoCall('cocoProblemGetBestValue', problem);