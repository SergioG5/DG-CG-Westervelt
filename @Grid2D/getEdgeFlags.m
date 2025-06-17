function edgePtrs = getEdgeFlags(G)
	numEdges = G.getNumEdges;
	nBndryEdges = G.getNumBoundaryEdges;
	edgePtrsFlag = false(numEdges, 1);
	for i = 1: nBndryEdges
		eList = G.bndry_.eList_(i);
		edgePtrsFlag(eList) = true;
	end
	nBndryEdges = nnz(edgePtrsFlag);
	nInteriorEdges = numEdges - nBndryEdges;
	edgePtrs = zeros(numEdges, 1);
	edgePtrs(edgePtrsFlag) = -(1: nBndryEdges);
	edgePtrs(~edgePtrsFlag) = 1: nInteriorEdges;
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