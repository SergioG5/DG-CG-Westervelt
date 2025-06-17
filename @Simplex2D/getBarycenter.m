function [G] = getBarycenter(T)
% ------------------------------------------------------------------------------
% G = getBarycenter() gets the barycenter coordinates of the simplex.
% ------------------------------------------------------------------------------
	G = transpose((T.v_(:,1) + T.v_(:,2) + T.v_(:,3))/3.0);
end

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