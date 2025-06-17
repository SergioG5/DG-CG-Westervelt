function [id] = getEdgeID(C, iSide)
% ----------------------------------------------------------------------------
% Returns the global edge id of side iSide. 
%
% ----------------------------------------------------------------------------
	id = C.eList_(iSide);
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