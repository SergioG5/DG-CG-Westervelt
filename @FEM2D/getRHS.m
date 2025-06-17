function [bf] = getRHS(fem, Thx, Tht, pde, nStep, uhN_1)
	p = fem.deg_;
	dim1DPp = p + 1;
	dim2DPp_3 = (p - 1)*(p - 2)/2;
	nInteriorNodes = Thx.getNumInteriorNodes();
	edgesPtrs = Thx.getEdgePtrs > 0;
	nInteriorEdges = nnz(edgesPtrs);
	numCells = Thx.getNumCells;
	N = nInteriorNodes + nInteriorEdges*(p-1) + numCells*dim2DPp_3;
	bf = zeros(dim1DPp*N, 1);
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	QW1D_ref = fem.op2D_.QW1D_;
	QX1D_ref = fem.op2D_.QX1D_;
	QW2D_ref = fem.op2D_.QW2D_;
	QX2D_ref = fem.op2D_.QX2D_;
	PXTB = fem.op2D_.PXTB_;
	PXTT = fem.op2D_.PXTT_;
	PXT = fem.op2D_.PXT_;
	QW_ref = fem.op2D_.QW_;
	Tht.getGeometry(nStep, Tt);
	QT = Tt.mapPoints(QX1D_ref); 
	ONES1D = ones(length(QW1D_ref), 1);
	ONES2D = ones(length(QW2D_ref), 1);
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1 : numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		DOF = fem.getDOF(Thx, kCell);
										% ------- pde source term  --------
		QX = Tx.mapPoints(QX2D_ref);
		QXT = [kron(ONES1D, QX), kron(QT, ONES2D)];
		FX = pde.f_(QXT(:, 1), QXT(:, 2), QXT(:, 3));
		bLocal = (Tt.getJacobian*Tx.getDetJac)*(PXT*(QW_ref.*FX));
		bf(DOF.ptrsInteriorDofsReduced) = bf(DOF.ptrsInteriorDofsReduced) ...
											+ bLocal(~DOF.bndryFlag);
		if(nStep == 1)
			U0X = pde.u0_(QX(:, 1), QX(:, 2));
			bottomLocal = Tx.getDetJac*(PXTB*(QW2D_ref.*U0X));
		else
			UHN_1X = PXTT'*uhN_1(DOF.ptrsGlobalDofs);
			bottomLocal = Tx.getDetJac*(PXTB*(QW2D_ref.*UHN_1X));
		end
		bf(DOF.ptrsInteriorDofsReduced) = bf(DOF.ptrsInteriorDofsReduced)...
										+ bottomLocal(~DOF.bndryFlag);
	end % end of loop over the elements
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