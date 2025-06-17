function uhN = computeL2proj(fem, Thx, Tht, F, nStep)
% ------------------------------------------------------------------------
% This function computes the L2 projection in the tensor product space of
% the function F
% ------------------------------------------------------------------------
	p = fem.degX_;
    q = fem.degT_;
	dimPp_3 = (p - 1)*(p - 2)/2;
	dim1DPq = q + 1;
	nInteriorNodes = Thx.getNumInteriorNodes();
	edgesPtrs = Thx.getEdgePtrs > 0;
	nInteriorEdges = nnz(edgesPtrs);
	numCells = Thx.getNumCells;
	nGlobalInteriorDofsX = nInteriorNodes + nInteriorEdges*(p - 1) ...
							+ numCells*dimPp_3;
	MAXNZ = 9 * numCells;
	M = spalloc(dim1DPq*nGlobalInteriorDofsX, ...
                    dim1DPq*nGlobalInteriorDofsX, MAXNZ);
	b = zeros(dim1DPq*nGlobalInteriorDofsX, 1);
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	Tht.getGeometry(nStep, Tt);
    % --------------------------------------------------------------------
    %                       LOAD OPERATORS
    % --------------------------------------------------------------------
    QT_ref = fem.op2D_.QX1D_;
    QX_ref = fem.op2D_.QX2D_;
    QW3D_ref = fem.op2D_.QW3D_;
    np1D = size(QT_ref, 1);
    np2D = size(QX_ref, 1);
    ONES1D = ones(np1D, 1);
    ONES2D = ones(np2D, 1);
    QT = Tt.mapPoints(QT_ref);
	hT = Tt.getJacobian;
    PTq_PX = fem.op2D_.PTQ_PX_;
    MT_MX = fem.op2D_.MT_MX_;
	% --------------------------------------------------------------------
	%					LOOP OVER THE MESH ELEMENTS
	% --------------------------------------------------------------------
    for iCell = 1: numCells
		Thx.getSimplex(iCell, Tx);
        detJ = Tx.getDetJac;
		Thx.getCell(iCell, cell);
		QX = Tx.mapPoints(fem.op2D_.QX2D_);
		QXT = [kron(ONES1D, QX), kron(QT, ONES2D)];
		DOF = fem.getDOF(Thx, iCell);
		Mlocal = (hT*detJ)*MT_MX;
		M(DOF.ptrsInteriorDofsReducedq, DOF.ptrsInteriorDofsReducedq) ...
			= M(DOF.ptrsInteriorDofsReducedq, DOF.ptrsInteriorDofsReducedq)  ...
			    + Mlocal(~DOF.bndryFlagq, ~DOF.bndryFlagq);
		FX = F(QXT(:, 1), QXT(:, 2), QXT(:, 3));
		bLocal = (hT*detJ)*(PTq_PX*(QW3D_ref.*FX));
		b(DOF.ptrsInteriorDofsReducedq) = b(DOF.ptrsInteriorDofsReducedq) ...
											+ bLocal(~DOF.bndryFlagq);
    end
    uhN = M\b;
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