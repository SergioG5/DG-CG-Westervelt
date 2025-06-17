function [bndryTag] = getBoundaryTag(C, iSide)
% ----------------------------------------------------------------------------
% Returns the boundary tag of boundary (edge) side iSide. 
%
% Warning: To avoid erroneous computations or possible misuse, the method 
%          aborts if that side is not a boundary face.
% ----------------------------------------------------------------------------

	isBoundary= (C.sList_(iSide) == 0);
	if isBoundary
		bndryTag = C.nList_(iSide);
	else
		error('Hmmm! asking boundary tag for an internal edge ...');
	end
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