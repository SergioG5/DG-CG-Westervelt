function [errGrad] = getErrorGradLinfDelta(fem, Thx, Tht, uhDelta1, uhDelta2, nStep)
	numCells = Thx.getNumCells;
	Tt    = Simplex1D();
	Tx    = Simplex2D();
	cell = Cell2D();
	errGrad = zeros(100, 1);
	PTLinfQ_DX = fem.op2D_.PTLinfQ_DX_;
	PTLinfQ_DY = fem.op2D_.PTLinfQ_DY_;
	QW2D_ref = fem.op2D_.QW2D_;
	Tht.getGeometry(nStep, Tt);
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1: numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		Ji = inv(Tx.getJacobian)';
		DOF = fem.getDOF(Thx, kCell);
		for s = 1: 100
			JDX = Ji(1, 1)*PTLinfQ_DX{s} + Ji(1, 2)*PTLinfQ_DY{s};
			JDY = Ji(2, 1)*PTLinfQ_DX{s} + Ji(2, 2)*PTLinfQ_DY{s};
			DUhXDelta1 = JDX'*uhDelta1(DOF.ptrsGlobalDofs);
			DUhXDelta2 = JDX'*uhDelta2(DOF.ptrsGlobalDofs);
			DUhYDelta1 = JDY'*uhDelta1(DOF.ptrsGlobalDofs);
			DUhYDelta2 = JDY'*uhDelta2(DOF.ptrsGlobalDofs);
			errGrad(s) = errGrad(s) ...
					+ (Tx.getDetJac)*dot(QW2D_ref, (DUhXDelta1 - DUhXDelta2).^2 ...
												 + (DUhYDelta1 - DUhYDelta2).^2);
		end
	end % end of loop over the elements
	errGrad = sqrt(max(errGrad));
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