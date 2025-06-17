function nodePtrs = getNodeFlags(G)
	numNodes = G.getNumNodes;
	nodePtrsFlag = false(numNodes, 1);
	nBndryEdges = length(G.bndry_.eList_);
	for i = 1: nBndryEdges
		eList = G.bndry_.eList_(i);
		eNodes = G.edges_.vList_(:, eList);
		nodePtrsFlag(eNodes) = true(2, 1);
	end
	nBndryNodes = nnz(nodePtrsFlag);
	nInteriorNodes = numNodes - nBndryNodes;
	nodePtrs = zeros(numNodes, 1);
	nodePtrs(nodePtrsFlag) = -(1: nBndryNodes);
	nodePtrs(~nodePtrsFlag) = 1: nInteriorNodes;
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