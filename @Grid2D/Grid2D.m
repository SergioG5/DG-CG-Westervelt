classdef Grid2D < handle

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Access = public )
	file_;  % file containing grid information
	nodes_; % structure containing the list of nodes
	nodePtrs_; %
	edgePtrs_; %
	edgeFlag_; % 
	edges_; % structure containing the list of edges
	cells_; % structure containing the list of cells
	bndry_; % 4 x numBndryEdges (uint32): list of boundary edges
	dg_;    % structure containing dg connectivity
end


% -----------------------------------------------------------------------------
%                         STATIC METHODS SECTION
% -----------------------------------------------------------------------------

methods (Static)

	[tag] = getDirichletTag();
	[tag] = getNeumannTag();
	[tag] = getRobinTag();

end

% -----------------------------------------------------------------------------
%                          METHODS SECTION
% -----------------------------------------------------------------------------
methods

function G = Grid2D(gridFile)
% ----------------------------------------------------------------------------
% Constructor:
%
% Grid2D(gridFile) gridFile is a binary file with grid information.
% 
% ----------------------------------------------------------------------------

	% ----------------------- reading grid file ------------------------------
	G.file_ = gridFile;
	fid = fopen(gridFile,'r');
	fprintf('Reading grid file "%s" ...\n',gridFile);
	% -------------------- READ VERTEX GEOMETRY ------------------------------
	G.nodes_.numNodes_ = fread(fid,[1 1],'*double');
	G.nodes_.xList_  = fread(fid,[2 G.nodes_.numNodes_],'*double');
	% -------------------- READ CELL CONNECTIVITY ----------------------------
	G.cells_.numCells_ = fread(fid,[1 1],'*double');
	G.cells_.vList_ = fread(fid,[3 G.cells_.numCells_],'*uint32');
	G.cells_.eList_ = fread(fid,[3 G.cells_.numCells_],'*uint32');
	G.cells_.mList_ = fread(fid,[1 G.cells_.numCells_],'*uint32');
	% -------------------- READ EDGE CONNECTIVITY ----------------------------
	G.edges_.numEdges_ = fread(fid,[1 1],'*double');
	G.edges_.vList_ = fread(fid,[2 G.edges_.numEdges_],'*uint32');
	G.edges_.nList_ = fread(fid,[2 G.edges_.numEdges_],'*uint32');
	G.edges_.sList_ = fread(fid,[2 G.edges_.numEdges_],'*uint8');
	G.edges_.hdgID_ = fread(fid,[1 G.edges_.numEdges_],'*uint32');
	G.edges_.length_ = zeros(1,G.edges_.numEdges_);
	G.edges_.numHdgEdges_ = max(G.edges_.hdgID_);
	G.edges_.globalHdgID_ = zeros(1,G.edges_.numHdgEdges_,'uint32');

	% ------------------- READ BOUNDARY EDGE CONNECTIVITY --------------------
	G.bndry_.numEdges_ = fread(fid,[1 1],'*double');
	G.bndry_.eList_ = [];
	G.bndry_.bType_ = [];
	if G.bndry_.numEdges_ > 0
		G.bndry_.eList_ = fread(fid,[1 G.bndry_.numEdges_],'*uint32');
		numTags = fread(fid,[1 1], '*double');
		G.bndry_.bType_ = fread(fid,[1 numTags],'*uint8');
	end
	fclose(fid);
	% ------------------------------------------------------------------------
	%          PROCESS SOME ADDITIONAL CONNECTIVITY INFORMATION
	% ------------------------------------------------------------------------
	G.nodePtrs_ = getNodeFlags(G);
	G.edgePtrs_ = getEdgeFlags(G);
	G.edgeFlag_ = G.edgePtrs_ > 0;
	disp('Processing additional connectivity information ... ')
	G.dg_.eList_ = zeros(3,G.cells_.numCells_,'uint32');
	G.dg_.nList_ = zeros(3,G.cells_.numCells_,'uint32');
	G.dg_.sList_ = zeros(3,G.cells_.numCells_,'uint8');
	G.dg_.rList_ = zeros(3,G.cells_.numCells_,'uint8');
	for iCell = 1:G.cells_.numCells_
		for iSide = 1:3
			iEdge = G.cells_.eList_(iSide,iCell);
			G.dg_.eList_(iSide,iCell) = G.edges_.hdgID_(iEdge);
			if G.edges_.nList_(1,iEdge) == iCell
				G.dg_.nList_(iSide,iCell) = G.edges_.nList_(2,iEdge);
				G.dg_.sList_(iSide,iCell) = G.edges_.sList_(2,iEdge);
				G.dg_.rList_(iSide,iCell) = 0; % POSITIVE ORIENTATION
			else
				G.dg_.nList_(iSide,iCell) = G.edges_.nList_(1,iEdge);
				G.dg_.sList_(iSide,iCell) = G.edges_.sList_(1,iEdge);
				G.dg_.rList_(iSide,iCell) = 1; % OPPOSITE ORIENTATION
			end
		end
	end
	numEdges = G.edges_.numEdges_;
	for k = 1:numEdges,
		i1 = G.edges_.vList_(1,k);
		i2 = G.edges_.vList_(2,k);
		G.edges_.length_(k) = norm(G.nodes_.xList_(:,i2) - ...
								G.nodes_.xList_(:,i1));
    end
end


	[numNodes] = getNumNodes(G);
	[node] = getNode(G, kNode);
	[numCells] = getNumCells(G);
	[materialId] = getMaterialId(G, kCell);
	[numEdges] = getNumEdges(G);
	[numEdges] = getNumInteriorEdges(G);
	[numEdges] = getNumBoundaryEdges(G);
	[] = getSimplex(G, kCell, T);
	[] = getCell(G, kCell, cell);
	[numInteriorNodes] = getNumInteriorNodes(G);
	[numBndryNodes] = getNumBndryNodes(G);
	[nodesPtrs] = getNodePtrs(G);
	[edgePtrs] = getEdgePtrs(G);
	nodePtrs = getNodeFlags(G);
	edgePtrs = getEdgeFlags(G);
	[hmin, hmax] = getGridSize(G);
	[QM] = evalQualityMeasures(G);
	[] = plot(G,varargin);

end					% END OF METHOD SECTION

% -----------------------------------------------------------------------------
%                          END OF CLASS DEFINITION
% -----------------------------------------------------------------------------
end

% ------------------------------------------------------------------------------
%                               END OF FILE
% ------------------------------------------------------------------------------








% -----------------------------------------------------------------------------
%                             Created by S. Gomez
%                          Department of Mathematics
%                         University of Pavia, Italia
%                                   (2020)
% -----------------------------------------------------------------------------