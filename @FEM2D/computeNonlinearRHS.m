function [bUh] = computeNonlinearRHS(fem, Thx, Tht, pde, nStep, uhNSmo)
	p = fem.degX_;
    q = fem.degT_;
	dimPp_3 = (p - 1)*(p - 2)/2;
	dim1DPq = q + 1;
	nInteriorNodes = Thx.getNumInteriorNodes();
	edgesPtrs = Thx.getEdgePtrs > 0;
	nInteriorEdges = nnz(edgesPtrs);
	numCells = Thx.getNumCells;
	nGlobalInteriorDofsX = nInteriorNodes + nInteriorEdges*(p - 1) ...
							+ numCells*dimPp_3;
	bUh = zeros(dim1DPq*nGlobalInteriorDofsX, 1);
	Tx = Simplex2D();
	Tt = Simplex1D();
	cell = Cell2D();
	Tht.getGeometry(nStep, Tt);

    % --------------------------------------------------------------------
    %                       LOAD OPERATORS
    % --------------------------------------------------------------------
	PXT = fem.op2D_.PTQ_PX_';
	DPXT = fem.op2D_.DPTQ_PX_';
	D2PXT = fem.op2D_.D2TQ_PX_';
    QW2D_ref = fem.op2D_.QW2D_;
    QW3D_ref = fem.op2D_.QW3D_;
	hT = Tt.getJacobian;
    PTQmo_PX = fem.op2D_.PTQmo_PX_;
    PTBottomQ_PX = fem.op2D_.PTBottomQ_PX_;
    DPTBottomQ_PX = fem.op2D_.DPTBottomQ_PX_;
	% -------------------------------------------------------------------
	%						LOOP OVER THE ELEMENTS
	% -------------------------------------------------------------------
	for kCell = 1: numCells
		Thx.getSimplex(kCell, Tx);
		Thx.getCell(kCell, cell);
		DOF = fem.getDOF(Thx, kCell);
        detJ = Tx.getDetJac;
		UHX = PXT*uhNSmo(DOF.ptrsGlobalDofs);
		DTUHX = DPXT*uhNSmo(DOF.ptrsGlobalDofs)/hT;
		D2TUHX = D2PXT*uhNSmo(DOF.ptrsGlobalDofs)/hT^2;
		% ----------------- CONTRIBUTION OF THE SOURCE TERM -------------
		bLocal = (hT*detJ)*(PTQmo_PX*(QW3D_ref.*(-pde.k_*DTUHX.^2 ...
												- pde.k_*UHX.*D2TUHX)));
		bUh(DOF.ptrsInteriorDofsReducedqmo) = ...
									bUh(DOF.ptrsInteriorDofsReducedqmo) ...
									+ bLocal(~DOF.bndryFlagqmo);

		% ---------------------------------------------------------------
		%	CONTRIBUTION OF THE INITIAL CONDITION OR PREVIOUS TIME STEP
		% ---------------------------------------------------------------
		UHNSmoX = PTBottomQ_PX'*uhNSmo(DOF.ptrsGlobalDofs);
		DTUHSmoX = DPTBottomQ_PX'*uhNSmo(DOF.ptrsGlobalDofs)/hT;
		bottomLocalWeak = -pde.k_*detJ*(PTBottomQ_PX*(QW2D_ref.*(UHNSmoX.*DTUHSmoX)));
        bUh(DOF.ptrsInteriorDofsReducedqmo) = ...
                bUh(DOF.ptrsInteriorDofsReducedqmo) ...
				+ bottomLocalWeak(~DOF.bndryFlagqmo);
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