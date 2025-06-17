function [isBndry] = isBoundaryEdge(C, iSide)
% ----------------------------------------------------------------------------
% Checks if the tetrahedron side iSide is a boundary face. 
%
% ----------------------------------------------------------------------------
	isBndry = (C.sList_(iSide) == 0);
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