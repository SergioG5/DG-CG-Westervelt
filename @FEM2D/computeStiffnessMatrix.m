function K = computeStiffnessMatrix(fem, Thx, Tht, pde, nStep)
% ------------------------------------------------------------------------
% [K] = computeStiffnessMatrix(Th, pde) computes the stiffness operator 
% in the 2D grid Th with the boundary conditions in the pde structure.
% ------------------------------------------------------------------------
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
	MAXNZ = 9 * numCells;
	K = spalloc(dim1DPq*nGlobalInteriorDofsX, ...
                    dim1DPq*nGlobalInteriorDofsX, MAXNZ);
	Tx = Simplex2D();
	Tt = Simplex1D();
	Cell = Cell2D();
	Tht.getGeometry(nStep, Tt);
	
    % --------------------------------------------------------------------
    %                       LOAD OPERATORS
    % --------------------------------------------------------------------
    QW2D_ref = fem.op2D_.QW2D_;
    QW3D_ref = fem.op2D_.QW3D_;
	hT = Tt.getJacobian;
    MXref = fem.op2D_.MX_;
	DX = fem.op2D_.DX_;
	DY = fem.op2D_.DY_;
    PTQ_DX = fem.op2D_.PTQ_DX_;
    PTQ_DY = fem.op2D_.PTQ_DY_;
    PTQmo_DX = fem.op2D_.PTQmo_DX_;
	DTQ_DX = fem.op2D_.DTQ_DX_;
	DTQ_DY = fem.op2D_.DTQ_DY_;
    PTQmo_DY = fem.op2D_.PTQmo_DY_;
    PTQmo_PX = fem.op2D_.PTQmo_PX_;
    D2Tq_PX = fem.op2D_.D2TQ_PX_;
    PTBqmo_PX = fem.op2D_.PTBottomQmo_PX_;
    DPTBq_PX = fem.op2D_.DPTBottomQ_PX_;

	% --------------------- LOOP OVER THE ELEMENTS ---------------------- %
	for iCell = 1: numCells
		Thx.getSimplex(iCell, Tx);
        detJ = Tx.getDetJac;
		Thx.getCell(iCell, Cell);
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
		D2Tlocal = (detJ/hT)*PTQmo_PX*(QW3D_ref.*D2Tq_PX') ...
                    + (detJ/hT)*PTBqmo_PX*(QW2D_ref.*DPTBq_PX');
		DampingLocal = detJ*pde.delta_*( ...
						  JDXqmo*(QW3D_ref.*JDTXq') ...
						+ JDYqmo*(QW3D_ref.*JDTYq') ...
					   );
		K(DOF.ptrsInteriorDofsReducedqmo, DOF.ptrsInteriorDofsReducedq) ...
			= K(DOF.ptrsInteriorDofsReducedqmo, DOF.ptrsInteriorDofsReducedq)  ...
			    + Klocal(~DOF.bndryFlagqmo, ~DOF.bndryFlagq) ...
			    + D2Tlocal(~DOF.bndryFlagqmo, ~DOF.bndryFlagq) ...
				+ DampingLocal(~DOF.bndryFlagqmo, ~DOF.bndryFlagq);

		if(nStep == 1)
			switch(fem.InitialCondition)
				case 'L2'
					KbottomLocal = detJ*MXref;
				case 'Ritz'
					JDX = Ji(1, 1)*DX + Ji(1, 2)*DY;
					JDY = Ji(2, 1)*DX + Ji(2, 2)*DY;
					KbottomLocal = detJ*( ...
						  JDX*(QW2D_ref.*JDX') ...
						+ JDY*(QW2D_ref.*JDY') ...
						);
				case 'Interpolation'
					KbottomLocal = speye(dim2DPp);
				otherwise
					error('Option not defined');
			end
		else
			KbottomLocal = detJ*MXref;
		end
        bottomDOFS = DOF.ptrsInteriorDofsReducedX ...
                        + dim1DPqmo*nGlobalInteriorDofsX;
        K(bottomDOFS, DOF.ptrsInteriorDofsReducedX) ...
	         = K(bottomDOFS, DOF.ptrsInteriorDofsReducedX) ...
              + KbottomLocal(~DOF.bndryFlagX, ~DOF.bndryFlagX);
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