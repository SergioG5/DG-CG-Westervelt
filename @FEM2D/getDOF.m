function [DOF] = getDOF(fem, Th, kCell)
	p = fem.degX_;
    q = fem.degT_;
	dim2DPp_3 = (p - 1)*(p - 2)/2;
	dim1DPq = q + 1;
    dim1DPqmo = q;
	dim2DPp = (p + 1)*(p + 2)/2;
	nLocalBubbleDofs = dim2DPp_3;
	nGlobalInteriorNodes = Th.getNumInteriorNodes;
	nGlobalNodes = Th.getNumNodes;
	edgesFlags = Th.edgeFlag_;
	nGlobalInteriorEdges = nnz(edgesFlags);
	nGlobalEdges = Th.getNumEdges;
	nCells = Th.getNumCells;
	nGlobalInteriorDofs = nGlobalInteriorNodes + nGlobalInteriorEdges*(p - 1) ...
							+ nCells*dim2DPp_3;
	nGlobalDofs = nGlobalNodes + nGlobalEdges*(p - 1) + nCells*dim2DPp_3;

	nodePtrs = Th.getNodePtrs();
	edgePtrs = Th.getEdgePtrs();
	kVerticesIndex = Th.cells_.vList_(:, kCell);
	kEdgesIndex = Th.cells_.eList_(:, kCell);
	bndryVerticesFlag = nodePtrs(kVerticesIndex) < 0;
	
	% ------------------------------------------------------------------------
	% 					LOCAL BOUNDARY DEGREES OF FREEDOM
	% ------------------------------------------------------------------------	
	DOF.bndryFlag = false(dim2DPp, 1);
	DOF.bndryFlag(1: 3) = bndryVerticesFlag;
	bEdgeFlags = edgePtrs(kEdgesIndex) < 0;
	cont = 4 + nLocalBubbleDofs;
	for mEdge = 1 : 3
		if(bEdgeFlags(mEdge))
			DOF.bndryFlag(cont: cont + p - 2) = true(p - 1, 1);
		end
		cont = cont + p - 1;
	end
	nLocalInteriorVertices = nnz(~bndryVerticesFlag);
	nLocalBndryVertices = nnz(bndryVerticesFlag);
	nLocalInteriorDofs = nnz(~DOF.bndryFlag);
	nLocalBndryDofs = nnz(DOF.bndryFlag);
	DOF.ptrsGlobalDofs = zeros(dim2DPp, 1);
	DOF.ptrsInteriorDofsReduced = zeros(nLocalInteriorDofs, 1);
	DOF.ptrsInteriorDofsFull = zeros(nLocalInteriorDofs, 1);
	DOF.ptrsBndryDofsFull = zeros(nLocalBndryDofs, 1);

	% ------------------------------------------------------------------------
	%						VERTEX DEGREES OF FREEDOM
	% ------------------------------------------------------------------------	
	DOF.ptrsGlobalDofs(1:3) = kVerticesIndex;
	DOF.ptrsInteriorDofsFull(1:nLocalInteriorVertices) = ...
											kVerticesIndex(~bndryVerticesFlag);
	DOF.ptrsInteriorDofsReduced(1:nLocalInteriorVertices) = ...
					nodePtrs(DOF.ptrsInteriorDofsFull(1:nLocalInteriorVertices));
	DOF.ptrsBndryDofsFull(1:nLocalBndryVertices) = ...
											kVerticesIndex(bndryVerticesFlag);
	% ------------------------------------------------------------------------
	% 						INTERIOR DEGREES OF FREEDOM
	% ------------------------------------------------------------------------	
	kInsideLimits = ((kCell-1)*nLocalBubbleDofs + 1: ...
					kCell*nLocalBubbleDofs) + nGlobalInteriorNodes;
	kInsideLimitsG = ((kCell-1)*nLocalBubbleDofs + 1: kCell*nLocalBubbleDofs) + nGlobalNodes;
	gLimits = nLocalInteriorVertices + 1: ...
							nLocalInteriorVertices +  length(kInsideLimits);
	DOF.ptrsInteriorDofsReduced(gLimits) = kInsideLimits;
	DOF.ptrsGlobalDofs(4: 3 + nLocalBubbleDofs) = kInsideLimitsG;
	DOF.ptrsInteriorDofsFull(gLimits) = kInsideLimitsG;

	% ------------------------------------------------------------------------
	% 						EDGE DEGREES OF FREEDOM
	% ------------------------------------------------------------------------	
	cont = nGlobalInteriorNodes + nCells*nLocalBubbleDofs;
	cont2 = nGlobalNodes + nCells*nLocalBubbleDofs;
	iLocalLimits = (nLocalInteriorVertices + length(kInsideLimitsG) + 1)...
        : (nLocalInteriorVertices + length(kInsideLimitsG) + p - 1);
	bLocalLimits = (nLocalBndryVertices + 1) : (nLocalBndryVertices + p - 1);
	mL = 3 + nLocalBubbleDofs;
	mLimits = mL + 1: mL + p - 1;
	for mEdge = 1: 3
		edge = kEdgesIndex(mEdge);
		orientation = Th.dg_.rList_(mEdge, kCell);
		if(~orientation)
			gEdge = edgePtrs(edge);
			eInteriorLimits = ((p-1)*(gEdge - 1) + 1: (p-1)*gEdge) + cont;
			eGlobalLimits = ((p-1)*(edge - 1) + 1: (p-1)*edge) + cont2;
		else
			gEdge = edgePtrs(edge);
			eInteriorLimits = ((p-1)*gEdge : -1 : (p-1)*(gEdge - 1) + 1) + cont;
			eGlobalLimits = ((p-1)*edge : -1 : (p-1)*(edge - 1) + 1) + cont2;
		end
		if(~bEdgeFlags(mEdge))
			DOF.ptrsInteriorDofsReduced(iLocalLimits) = eInteriorLimits;
			DOF.ptrsInteriorDofsFull(iLocalLimits) = eGlobalLimits;
			iLocalLimits = iLocalLimits + p - 1;
		else
			DOF.ptrsBndryDofsFull(bLocalLimits) = ...
							eGlobalLimits;
			bLocalLimits = bLocalLimits + p - 1;
		end
		DOF.ptrsGlobalDofs(mLimits) = eGlobalLimits;
		mLimits = mLimits + p - 1;
	end
	% ------------------------------------------------------------------------
	% 						SPACE-TIME TENSOR PRODUCT
	% ------------------------------------------------------------------------	
	aux.bndryFlagq = false(dim1DPq*dim2DPp, 1);
    aux.bndryFlagqmo = false(dim1DPqmo*dim2DPp, 1);
	aux.ptrsGlobalDofs = zeros(dim1DPq*dim2DPp, 1);
	aux.ptrsInteriorDofsFull = zeros(dim1DPq*nLocalInteriorDofs, 1);
	aux.ptrsInteriorDofsReducedq = zeros(dim1DPq*nLocalInteriorDofs, 1);
    aux.ptrsInteriorDofsReducedqmo = zeros(dim1DPqmo*nLocalInteriorDofs, 1);
	aux.ptrsBndryDofsReduced = zeros(nnz(DOF.bndryFlag), 1);
    aux.ptrsInteriorDofsReducedX = DOF.ptrsInteriorDofsReduced;
    aux.bndryFlagX = DOF.bndryFlag;
    for j = 1: q + 1
		jLimits = (j - 1)*dim2DPp + 1: j*dim2DPp;
		jInteriorLimits = (j - 1)*nLocalInteriorDofs + 1: j*nLocalInteriorDofs;
		jBndryLimits = (j - 1)*nLocalBndryDofs + 1: j*nLocalBndryDofs;
		aux.bndryFlagq(jLimits) = DOF.bndryFlag;
		aux.ptrsGlobalDofs(jLimits) = DOF.ptrsGlobalDofs + (j - 1)*nGlobalDofs;
		aux.ptrsInteriorDofsFull(jInteriorLimits) = DOF.ptrsInteriorDofsFull...
											+ (j - 1)*nGlobalDofs;
		aux.ptrsInteriorDofsReducedq(jInteriorLimits) = DOF.ptrsInteriorDofsReduced ...
												+ (j - 1)*nGlobalInteriorDofs;
		aux.ptrsBndryDofsFull(jBndryLimits) = DOF.ptrsBndryDofsFull ...
											+ (j - 1)*nGlobalDofs;
    end
    for j = 1: q
        jLimits = (j - 1)*dim2DPp + 1: j*dim2DPp;
		jInteriorLimits = (j - 1)*nLocalInteriorDofs + 1: j*nLocalInteriorDofs;
		aux.ptrsInteriorDofsReducedqmo(jInteriorLimits) ...
            = DOF.ptrsInteriorDofsReduced + (j - 1)*nGlobalInteriorDofs;
        aux.bndryFlagqmo(jLimits) = DOF.bndryFlag;
    end
	DOF = aux;
end  % end function

% -----------------------------------------------------------------------------
% Created by 
%
% Sergio Gomez, sergio.gomezmacias@unimib.it
% Department of Mathematics and Applications
% University of Milano-Bicocca (UNIMIB)
%
%                                   (2025)
% -----------------------------------------------------------------------------