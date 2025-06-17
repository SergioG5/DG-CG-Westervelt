classdef Cell2D < handle
% ----------------------------------------------------------------------------
% Definition of a triangular cell connectivity for discontinuous finite element 
% computations. Vertices are oriented counter clockwise.
%
% Sides (local edges):
%
%             Side 1: (2,3) opposite to node 1
%             Side 2: (3,1) opposite to node 2
%             Side 3: (1,2) opposite to node 3
%
% Triangle side i has a pointer to a global edge and a local orientation 
% relative to that global edge. Each side is, potentially, shared by a 
% neighboring cell. That neighbor sees the current triangle in its local side.
% Otherwise that edge is a boundary edge.
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Access = ?Grid2D )
	eList_;		% 3 x 1 (uint32): edge IDs 
	nList_;		% 3 x 1 (uint32): neighbor ID
	sList_;		% 3 x 1 (uint8): neighbor local side ID
	oList_;		% 3 x 1 (uint8): edge orientation

end


% -----------------------------------------------------------------------------
%                          METHODS SECTION
% -----------------------------------------------------------------------------
methods


function [C] = Cell2D()
% ----------------------------------------------------------------------------
% Constructor of an object of class Cell2D
%
% ----------------------------------------------------------------------------
	C.eList_ = zeros(3,1,'uint32');
	C.nList_ = zeros(3,1,'uint32');
	C.sList_ = zeros(3,1,'uint8');
	C.oList_ = zeros(3,1,'uint8');
end

	[id] = getEdgeID(C, iSide);
	[id] = getNeighborID(C, iSide);
	[id] = getSideID(C, iSide);
	[id] = getOrientation(C, iSide);
	[isBndry] = isBoundaryEdge(C, iSide);
	[bndryTag] = getBoundaryTag(C, iSide);

end					% END OF METHOD SECTION

% -----------------------------------------------------------------------------
%                          END OF CLASS DEFINITION
% -----------------------------------------------------------------------------
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