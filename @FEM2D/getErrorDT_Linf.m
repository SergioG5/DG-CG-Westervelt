function [errDT_Linf] = getErrorDT_Linf(fem, Thx, Tht, uh, Func, nStep)
	numCells = Thx.getNumCells;
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	QW2D_ref = fem.op2D_.QW2D_;
	QX2D_ref = fem.op2D_.QX2D_;
	Tht.getGeometry(nStep, Tt);
	tLimits = Tt.getVertices;
	tArray = linspace(tLimits(1), tLimits(2), 100);
	hT = Tt.getDiameter;
	DPTLinfQ_PX = fem.op2D_.DPTLinfQ_PX_;
	errDT_Linf = zeros(100, 1);
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1 : numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		DOF = fem.getDOF(Thx, kCell);
										% ------- pde source term  --------
		QX = Tx.mapPoints(QX2D_ref);
		for s = 1: 100
			DT_UX = Func(QX(:, 1), QX(:, 2), tArray(s));
			DT_UhX = (1/hT)*DPTLinfQ_PX{s}'*uh(DOF.ptrsGlobalDofs);
			errDT_Linf(s) = errDT_Linf(s) + (Tx.getDetJac)*dot(QW2D_ref, (DT_UX - DT_UhX).^2);
		end
	end % end of loop over the elements
	errDT_Linf = sqrt(max(errDT_Linf));
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