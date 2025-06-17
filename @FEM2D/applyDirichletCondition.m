function uhGlobal = applyDirichletCondition(fem, Thx, Tht, uh, pde, nStep)
% ----------------------------------------------------------------------
% This function set the boundary degrees of freedom of the solution uh
% obtained from the reduced system
% ----------------------------------------------------------------------
	% ----------------- AUXILIARY VARIABLES AND OPERATORS --------------
	p = fem.degX_;
    q = fem.degT_;
	dim2DPp_3 = (p - 1)*(p - 2)/2;
	dim2DPp = (p + 1)*(p + 2)/2;
	dim1DPq = q + 1;
	nNodes = Thx.getNumNodes();
	nEdges = Thx.getNumEdges();
	nCells = Thx.getNumCells;
	Tt = Simplex1D;
	Tx = Simplex2D;
	cell = Cell2D();
	N = nNodes + nEdges*(p - 1) + nCells*dim2DPp_3;
	lagrangePoints = fem.op2D_.BX_.getLagrangePoints;
	legendrePoints = fem.op2D_.BTq_.interpolationNodes;
	Tht.getGeometry(nStep, Tt);
	vT = Tt.getVertices;
	tNmo = vT(1);
	tN = vT(2);
	uhGlobal = zeros(dim1DPq*N, 1);
	% ----------------------------- QUADRATURES ------------------------
	QT_ref = fem.op2D_.QX1D_;
	QW1D_ref = fem.op2D_.QW1D_;
    np1D = size(QT_ref, 1);
    QT = Tt.mapPoints(QT_ref);
	hT = Tt.getJacobian;
	projMatrix = [
		fem.op2D_.PTBottomQ_';
		(1/hT)*fem.op2D_.DPTTopQ_';
		(fem.op2D_.PTQm2_*(QW1D_ref.*fem.op2D_.DTQ_'));
	];
	% -----------------------------------------------------------------
	%						CYCLE ON THE MESH ELEMENTS
	% -----------------------------------------------------------------
	for kCell = 1 : nCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		DOF = fem.getDOF(Thx, kCell);
		xNodes = Tx.mapPoints(lagrangePoints);
		% ----------------- DIFFERENT WAYS TO ENFORCE THE BC ---------
		switch(fem.BoundaryCondition)
			case 'Interpolation'	% LAGRANGE INTERPOLATION
				tNodes = Tt.mapPoints(legendrePoints);
				xtNodes = [kron(ones(dim1DPq, 1), xNodes), kron(tNodes, ...
							ones(dim2DPp, 1))];
				gD = pde.u_(xtNodes(:, 1), xtNodes(:, 2), xtNodes(:, 3));
			case 'Projection'		% WALKINGTON PROJECTION
				gDT = zeros(np1D, dim2DPp);
				for s = 1: dim2DPp
					gDT(:, s) = pde.ut_(xNodes(s, 1), xNodes(s, 2), QT);
				end
				bD = [
					pde.u_(xNodes(:, 1), xNodes(:, 2), tNmo)';
					pde.ut_(xNodes(:, 1), xNodes(:, 2), tN)';
					hT*(fem.op2D_.PTQm2_*(QW1D_ref.*gDT))
				];
				proj = projMatrix\bD;
				proj = proj';
				gD = proj(:);
		end
		uhGlobal(DOF.ptrsInteriorDofsFull) = ...
										uh(DOF.ptrsInteriorDofsReducedq);
		uhGlobal(DOF.ptrsBndryDofsFull) = gD(DOF.bndryFlagq);
	end
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