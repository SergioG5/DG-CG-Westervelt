function [errL2T] = getErrorL2T(fem, Thx, Tht, uh, Func, nStep)
	numCells = Thx.getNumCells;
	Tx = Simplex2D();
	Tt = Simplex1D();
	QW2D_ref = fem.op2D_.QW2D_;
	QX2D_ref = fem.op2D_.QX2D_;
	PXTTopQ = fem.op2D_.PTTopQ_PX_';
	Tht.getGeometry(nStep, Tt);
	tf = max(Tt.getVertices);
	errL2T = 0.0;
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1 : numCells
		Thx.getSimplex(kCell, Tx);
		DOF = fem.getDOF(Thx, kCell);
										% ------- pde source term  --------
		QX = Tx.mapPoints(QX2D_ref);
		UX = Func(QX(:, 1), QX(:, 2), tf);
		UhX = PXTTopQ*uh(DOF.ptrsGlobalDofs);
		errL2T = errL2T + (Tx.getDetJac)*dot(QW2D_ref, (UX - UhX).^2);
	end % end of loop over the elements
	errL2T = sqrt(errL2T);
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