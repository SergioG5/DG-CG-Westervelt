function [b] = computeRHS(fem, Thx, Tht, pde, nStep, uhN_1)
	p = fem.degX_;
    q = fem.degT_;
	dimPp_3 = (p - 1)*(p - 2)/2;
	dim2DPp = (p + 1)*(p + 2)/2;
	dim1DPq = q + 1;
    dim1DPqmo = q;
	nInteriorNodes = Thx.getNumInteriorNodes();
	edgesPtrs = Thx.getEdgePtrs > 0;
	nInteriorEdges = nnz(edgesPtrs);
	numCells = Thx.getNumCells;
	nGlobalInteriorDofsX = nInteriorNodes + nInteriorEdges*(p - 1) ...
							+ numCells*dimPp_3;
	lagrangePointsX = fem.op2D_.BX_.getLagrangePoints;
	lagrangePointsT = fem.op2D_.BTq_.interpolationNodes;
	b = zeros(dim1DPq*nGlobalInteriorDofsX, 1);
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	xDim = 1: dim2DPp;
	yDim = xDim + dim2DPp;
	Tht.getGeometry(nStep, Tt);
	vT = Tt.getVertices;
	tNmo = vT(1);
	tN = vT(2);

    % --------------------------------------------------------------------
    %                       LOAD OPERATORS
    % --------------------------------------------------------------------
    QT_ref = fem.op2D_.QX1D_;
	QW1D_ref = fem.op2D_.QW1D_;
    QW2D_ref = fem.op2D_.QW2D_;
    QX_ref = fem.op2D_.QX2D_;
    QW3D_ref = fem.op2D_.QW3D_;
    np1D = size(QT_ref, 1);
    np2D = size(QX_ref, 1);
    ONES1D = ones(np1D, 1);
    ONES2D = ones(np2D, 1);
    QT = Tt.mapPoints(QT_ref);
	hT = Tt.getJacobian;
    PX = fem.op2D_.PX_;
    MXref = fem.op2D_.MX_;
	DX = fem.op2D_.DX_;
	DY = fem.op2D_.DY_;
    PTQ_DX = fem.op2D_.PTQ_DX_;
    PTQ_DY = fem.op2D_.PTQ_DY_;
	DTQ_DX = fem.op2D_.DTQ_DX_;
	DTQ_DY = fem.op2D_.DTQ_DY_;
    PTQmo_DX = fem.op2D_.PTQmo_DX_;
    PTQmo_DY = fem.op2D_.PTQmo_DY_;
    PTQmo_PX = fem.op2D_.PTQmo_PX_;
    D2TQ_PX = fem.op2D_.D2TQ_PX_;
    PTTopQ_PX = fem.op2D_.PTTopQ_PX_;
    DPTTopQ_PX = fem.op2D_.DPTTopQ_PX_;
    PTBottomQmo_PX = fem.op2D_.PTBottomQmo_PX_;
    DPTBottomQ_PX = fem.op2D_.DPTBottomQ_PX_;
	% -------------------------------------------------------------------
	%						LOOP OVER THE ELEMENTS
	% -------------------------------------------------------------------
	for iCell = 1: numCells
		Thx.getSimplex(iCell, Tx);
        detJ = Tx.getDetJac;
		Thx.getCell(iCell, cell);
		QX = Tx.mapPoints(fem.op2D_.QX2D_);
		QXT = [kron(ONES1D, QX), kron(QT, ONES2D)];
		Ji = inv(Tx.getJacobian)';
		DOF = fem.getDOF(Thx, iCell);
		JDXq = Ji(1, 1)*PTQ_DX + Ji(1, 2)*PTQ_DY;
		JDYq = Ji(2, 1)*PTQ_DX + Ji(2, 2)*PTQ_DY;
		JDTXq = Ji(1, 1)*DTQ_DX + Ji(1, 2)*DTQ_DY;
		JDTYq = Ji(2, 1)*DTQ_DX + Ji(2, 2)*DTQ_DY;
        JDXqmo = Ji(1, 1)*PTQmo_DX + Ji(1, 2)*PTQmo_DY;
		JDYqmo = Ji(2, 1)*PTQmo_DX + Ji(2, 2)*PTQmo_DY;
		Klocal =   (hT*detJ*pde.c_^2)*( ...
					  JDXqmo*(QW3D_ref.*JDXq') ...
					+ JDYqmo*(QW3D_ref.*JDYq') ...
				 );
		D2Tlocal = (detJ/hT)*PTQmo_PX*(QW3D_ref.*D2TQ_PX') ...
                    + (detJ/hT)*PTBottomQmo_PX*(QW2D_ref.*DPTBottomQ_PX');
		DampingLocal = detJ*pde.delta_*( ...
						  JDXqmo*(QW3D_ref.*JDTXq') ...
						+ JDYqmo*(QW3D_ref.*JDTYq') ...
					   );
		xNodes = Tx.mapPoints(lagrangePointsX);
		tNodes = Tt.mapPoints(lagrangePointsT);
		% --------------- APPROXIMATION OF THE DIRICHLET DATA -----------
		switch(fem.BoundaryCondition)
			case 'Interpolation'
				xtNodes = [kron(ones(dim1DPq, 1), xNodes), ...
                    kron(tNodes, ones(dim2DPp, 1))];
				gD = pde.u_(xtNodes(:, 1), xtNodes(:, 2), xtNodes(:, 3));
			case 'Projection'
				projMatrix = [
					fem.op2D_.PTBottomQ_';
					(1/hT)*fem.op2D_.DPTTopQ_';
					(fem.op2D_.PTQm2_*(QW1D_ref.*fem.op2D_.DTQ_'));
				];
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
		b(DOF.ptrsInteriorDofsReducedqmo) = b(DOF.ptrsInteriorDofsReducedqmo) ...
			- Klocal(~DOF.bndryFlagqmo, DOF.bndryFlagq)*gD(DOF.bndryFlagq) ...
			- D2Tlocal(~DOF.bndryFlagqmo, DOF.bndryFlagq)*gD(DOF.bndryFlagq) ...
			- DampingLocal(~DOF.bndryFlagqmo, DOF.bndryFlagq)*gD(DOF.bndryFlagq);

		% ----------------- CONTRIBUTION OF THE SOURCE TERM -------------
		FX = pde.f_(QXT(:, 1), QXT(:, 2), QXT(:, 3));
		bLocal = (hT*detJ)*(PTQmo_PX*(QW3D_ref.*FX));
		b(DOF.ptrsInteriorDofsReducedqmo) = b(DOF.ptrsInteriorDofsReducedqmo) ...
											+ bLocal(~DOF.bndryFlagqmo);

		% ---------------------------------------------------------------
		%	CONTRIBUTION OF THE INITIAL CONDITION OR PREVIOUS TIME STEP
		% ---------------------------------------------------------------
		if(nStep == 1)
			U0X = pde.u0_(QX(:, 1), QX(:, 2));
			switch(fem.InitialCondition)
				case 'L2'
					KbottomLocal = detJ*MXref;
					bottomLocal = detJ*(PX*(QW2D_ref.*U0X));
				case 'Ritz'
					JDX = kron(Ji(:, 1), DX) + kron(Ji(:, 2), DY);
					DX_U0 = pde.u0x_(QX(:, 1), QX(:, 2));
					DY_U0 = pde.u0y_(QX(:, 1), QX(:, 2));
					KbottomLocal = detJ*( ...
						  JDX(xDim, :)*(QW2D_ref.*JDX(xDim, :)') ...
						+ JDX(yDim, :)*(QW2D_ref.*JDX(yDim, :)') ...
						);
					bottomLocal = detJ*(JDX(xDim, :)*(QW2D_ref.*DX_U0) ...
								+ JDX(yDim, :)*(QW2D_ref.*DY_U0));
				case 'Interpolation'
					KbottomLocal = speye(dim2DPp);
					bottomLocal = pde.u0_(xNodes(:, 1), xNodes(:, 2));
				otherwise
					error('Option not defined');
			end
            V0X = pde.v0_(QX(:, 1), QX(:, 2));
            bottomLocalWeak = detJ*(PTBottomQmo_PX*(QW2D_ref.*(1 + pde.k_*U0X).*V0X));
		else
			UHN_1X = PTTopQ_PX'*uhN_1(DOF.ptrsGlobalDofs);
			bottomLocal = detJ*(PX*(QW2D_ref.*UHN_1X));
            DTUHN_1X = DPTTopQ_PX'*uhN_1(DOF.ptrsGlobalDofs);
			bottomLocalWeak = (detJ/hT)*(PTBottomQmo_PX*(QW2D_ref.*(1 + pde.k_*UHN_1X).*DTUHN_1X));
			KbottomLocal = detJ*MXref;
		end
        b(DOF.ptrsInteriorDofsReducedqmo) = ...
                b(DOF.ptrsInteriorDofsReducedqmo) ...
				+ bottomLocalWeak(~DOF.bndryFlagqmo);

		% We assume that we interpolate at the initial time!
        bottomDOFS = DOF.ptrsInteriorDofsReducedX ...
                        + dim1DPqmo*nGlobalInteriorDofsX;
		gD = pde.u_(xNodes(:, 1), xNodes(:, 2), tNmo);
 		b(bottomDOFS) = b(bottomDOFS) + bottomLocal(~DOF.bndryFlagX) ...
		- KbottomLocal(~DOF.bndryFlagX, DOF.bndryFlagX)*gD(DOF.bndryFlagX);
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