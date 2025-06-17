function [errGrad] = getErrorGradLinf(fem, Thx, Tht, uh, ux, uy, nStep)
	p = fem.degX_;
	q = fem.degT_;
	dim2DPp = (p + 1)*(p + 2)/2;
	dim1DPq = q + 1;
	numCells = Thx.getNumCells;
	Tt    = Simplex1D();
	Tx    = Simplex2D();
	cell = Cell2D();
	errGrad = zeros(100, 1);
	PTLinfQ_DX = fem.op2D_.PTLinfQ_DX_;
	PTLinfQ_DY = fem.op2D_.PTLinfQ_DY_;
	QX2D_ref = fem.op2D_.QX2D_;
	QW2D_ref = fem.op2D_.QW2D_;
	Tht.getGeometry(nStep, Tt);
	tLimits = Tt.getVertices;
	tArray = linspace(tLimits(1), tLimits(2), 100);
	xDim = 1: dim1DPq*dim2DPp;
	yDim = xDim + dim1DPq*dim2DPp;
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1: numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		Ji = inv(Tx.getJacobian)';
		DOF = fem.getDOF(Thx, kCell);
		QX = Tx.mapPoints(QX2D_ref);
		for s = 1: 100
			JDX = Ji(1,1)*PTLinfQ_DX{s} + Ji(1,2)*PTLinfQ_DY{s};
			JDY = Ji(2,1)*PTLinfQ_DX{s} + Ji(2,2)*PTLinfQ_DY{s};
			DUX = ux(QX(:, 1), QX(:, 2), tArray(s));
			DUY = uy(QX(:, 1), QX(:, 2), tArray(s));
			DUhX = JDX'*uh(DOF.ptrsGlobalDofs);
			DUhY = JDY'*uh(DOF.ptrsGlobalDofs);
			errGrad(s) = errGrad(s) ...
					+ (Tx.getDetJac)*dot(QW2D_ref, (DUhX - DUX).^2 ...
												 + (DUhY - DUY).^2);
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