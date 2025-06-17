function [errGrad] = getErrorGrad(fem, Thx, Tht, uh, ux, uy, nStep)
	p = fem.degX_;
	q = fem.degT_;
	dim2DPp = (p + 1)*(p + 2)/2;
	dim1DPq = q + 1;
	numCells = Thx.getNumCells;
	Tt    = Simplex1D();
	Tx    = Simplex2D();
	cell = Cell2D();
	errGrad = 0.0;
	PTQ_DX = fem.op2D_.PTQ_DX_;
	PTQ_DY = fem.op2D_.PTQ_DY_;
	QW1D_ref = fem.op2D_.QW1D_;
	QX1D_ref = fem.op2D_.QX1D_;
	QW2D_ref = fem.op2D_.QW2D_;
	QX2D_ref = fem.op2D_.QX2D_;
	QW_ref = fem.op2D_.QW3D_;
	Tht.getGeometry(nStep, Tt);
	QT = Tt.mapPoints(QX1D_ref);
	ONES1D = ones(length(QW1D_ref), 1);
	ONES2D = ones(length(QW2D_ref), 1);
	xDim = 1: dim1DPq*dim2DPp;
	yDim = xDim + dim1DPq*dim2DPp;
	% -------------- LOOP OVER THE ELEMENTS -------------- %
	for kCell = 1: numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		Ji = inv(Tx.getJacobian)';
		JDX = kron(Ji(:,1), PTQ_DX) + kron(Ji(:,2), PTQ_DY);
		DOF = fem.getDOF(Thx, kCell);
		QX = Tx.mapPoints(QX2D_ref);
		QXT = [kron(ONES1D, QX), kron(QT, ONES2D)];
		DUX = ux(QXT(:, 1), QXT(:, 2), QXT(:, 3));
		DUY = uy(QXT(:, 1), QXT(:, 2), QXT(:, 3));
		DUhX = JDX(xDim, :)'*uh(DOF.ptrsGlobalDofs);
		DUhY = JDX(yDim, :)'*uh(DOF.ptrsGlobalDofs);
		errGrad = errGrad ...
				+ (Tt.getJacobian*Tx.getDetJac)*dot(QW_ref, (DUhX - DUX).^2 ...
												 + (DUhY - DUY).^2);
	end % end of loop over the elements
	errGrad = sqrt(errGrad);
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