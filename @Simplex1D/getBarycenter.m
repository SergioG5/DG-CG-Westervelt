function [barycenter] = getBarycenter(T)
% ------------------------------------------------------------------------------
% Get barycenter of 1D-simplex.
% OUTPUT:
%   midpoint of current cell.
% ------------------------------------------------------------------------------
	barycenter = (T.v_(1) + T.v_(2))/2.0;
end


% ------------------------------------------------------------------------------
%                               END OF FILE
% ------------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Created by 
%
% Paul Castillo, paul.castillo@upr.edu
% Department of Mathematical Sciences 
% University of Puerto Rico, Mayaguez Campus (UPRM)
%
% Sergio Gomez, sergio.gomezmacias@unimib.it
% Department of Mathematics and Applications
% University of Milano-Bicocca (UNIMIB)
%
%                                   (2020)
% -----------------------------------------------------------------------------