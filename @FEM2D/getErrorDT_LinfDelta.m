function [errDT_Linf] = getErrorDT_LinfDelta(fem, Thx, Tht, uhDelta1, uhDelta2, nStep)
	numCells = Thx.getNumCells;
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	QW2D_ref = fem.op2D_.QW2D_;
	Tht.getGeometry(nStep, Tt);
	hT = Tt.getDiameter;
	DPTLinfQ_PX = fem.op2D_.DPTLinfQ_PX_;
	errDT_Linf = zeros(100, 1);
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1 : numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		DOF = fem.getDOF(Thx, kCell);
										% ------- pde source term  --------
		for s = 1: 100
			DT_UhXDelta1 = (1/hT)*DPTLinfQ_PX{s}'*uhDelta1(DOF.ptrsGlobalDofs);
			DT_UhXDelta2 = (1/hT)*DPTLinfQ_PX{s}'*uhDelta2(DOF.ptrsGlobalDofs);
			errDT_Linf(s) = errDT_Linf(s) + (Tx.getDetJac)*dot(QW2D_ref, (DT_UhXDelta1 - DT_UhXDelta2).^2);
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