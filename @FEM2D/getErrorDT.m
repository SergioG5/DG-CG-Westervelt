function [errL2] = getErrorDT(fem, Thx, Tht, uh, Func, nStep)
	numCells = Thx.getNumCells;
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	QW1D_ref = fem.op2D_.QW1D_;
	QX1D_ref = fem.op2D_.QX1D_;
	QW2D_ref = fem.op2D_.QW2D_;
	QX2D_ref = fem.op2D_.QX2D_;
	QW_ref = fem.op2D_.QW3D_;
	DPTQ_PX = fem.op2D_.DPTQ_PX_';
	Tht.getGeometry(nStep, Tt);
	QT = Tt.mapPoints(QX1D_ref); 
	ONES1D = ones(length(QW1D_ref), 1);
	ONES2D = ones(length(QW2D_ref), 1);
    hT = Tt.getDiameter;
	errL2 = 0.0;
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1 : numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		DOF = fem.getDOF(Thx, kCell);
										% ------- pde source term  --------
		QX = Tx.mapPoints(QX2D_ref);
		QXT = [kron(ONES1D, QX), kron(QT, ONES2D)];
		UX = Func(QXT(:, 1), QXT(:, 2), QXT(:, 3));
		DTUhX = (1/hT)*DPTQ_PX*uh(DOF.ptrsGlobalDofs);
		errL2 = errL2 + (hT*Tx.getDetJac)*dot(QW_ref, (UX - DTUhX).^2);
	end % end of loop over the elements
	errL2 = sqrt(errL2);
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