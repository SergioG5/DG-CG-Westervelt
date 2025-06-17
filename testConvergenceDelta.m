clear all
close all
clc
% -------------------------------------------------------------
% 					TEST PARAMETERS 					
% -------------------------------------------------------------
testNumber = 1;
deltaArray = 10.^(-(2:2:10));
pMin       = 1;
pMax       = 3;
q = 4;
nMeshRefinements = length(deltaArray);
ICopt = 'Ritz';
BCopt = 'Projection';

tol = 1.0e-12;
nMaxIt = 15;

% --------------------- MESHES ----------------------------- %

GRID = {'grid-10x10.mat'};
N = 20;

% ------------------------------------------------------------- %
% 					ERROR VECTORS  					
% ------------------------------------------------------------- %
nDOFS = zeros(nMeshRefinements, pMax);
errL2 = zeros(nMeshRefinements, pMax);
errDT = zeros(nMeshRefinements, pMax);
errDT_L2T = zeros(nMeshRefinements, pMax);
errDT_Linf = zeros(nMeshRefinements, pMax);
errGrad_Linf = zeros(nMeshRefinements, pMax);
errL2T = zeros(nMeshRefinements, pMax);
errL2H1 = zeros(nMeshRefinements, pMax);

% ------------------------------------------------------------- %
% 					TESTS DESCRIPTION  					
% ------------------------------------------------------------- % 
switch(testNumber)
	case 1
		testName = "Gomez-Meliani";
        pde.c_ = 1.0;
		pde.k_ = 0.3;
		A = 1.0e-2;
		pde.u_ = @(x, y, t) zeros(size(x));
        pde.ut_ = @(x, y, t) zeros(size(x));
		pde.ux_ = @(x, y, t) zeros(size(x));
		pde.uy_ = @(x, y, t) zeros(size(x));
		pde.u0_ = @(x, y) A*sin(pi*x).*sin(pi*y);
		pde.u0x_ = @(x, y) A*pi*cos(pi*x).*sin(pi*y);
		pde.u0y_ = @(x, y) A*pi*cos(pi*y).*sin(pi*x);
        pde.v0_ = @(x, y) sin(pi*x).*sin(pi*y);
		pde.f_ = @(x, y, t) zeros(size(x));
		pde.tf_ = 1.0;
	case 2
		testName = "Steepening-wave-front";
        pde.c_ = 2000;
		pde.k_ = -10;
		a = 400;
		alpha = 5.0e4;
		sigma = 3.0e-2;
		pde.u_ = @(x, y, t) zeros(size(x));
        pde.ut_ = @(x, y, t) zeros(size(x));
		pde.ux_ = @(x, y, t) zeros(size(x));
		pde.uy_ = @(x, y, t) zeros(size(x));
		pde.u0_ = @(x, y) zeros(size(x));
		pde.u0x_ = @(x, y) zeros(size(x));
		pde.u0y_ = @(x, y) zeros(size(x));
        pde.v0_ = @(x, y) zeros(size(x));
		pde.f_ = @(x, y, t) (a/sqrt(sigma))*exp(-alpha*t).*exp(-((x - 0.5).^2 + (y - 0.5).^2)/(2*sigma^2));
		pde.tf_ = 2.0e-4;
end

% ------------------------------------------------------------- %
% 					MAIN LOOP  					
% ------------------------------------------------------------- %

% ------------------LOOP OVER DEGREES --------------------- %
for p = pMin : pMax
	dim = (p + 1)*(p + 2)/2;
    fem = FEM2D(p, q, ICopt, BCopt);
	Grid= strcat('Grids/Dirichlet/',GRID{1});
	ThX = Grid2D(Grid);
	tPoints = linspace(0, pde.tf_, N(1) + 1);
	ht = pde.tf_/N(1);
	ThT = Grid1D(tPoints');
	nSteps = N(1);
	% -------------- SOLUTION OF VANISHING LIMIT -------------------- %
	pde.delta_ = 0.0;
	uhGlobal0 = cell(nSteps, 1);
	for n = 1 : nSteps
		if(n == 1)
			K = fem.computeStiffnessMatrix(ThX, ThT, pde, n);
			b = fem.computeRHS(ThX, ThT, pde, n);
			uhSmo = K\b;
		else
			uhSmo = uh;
			if(n == 2)
				K = fem.computeStiffnessMatrix(ThX, ThT, pde, n);
				K = decomposition(K);
			end
			b = fem.computeRHS(ThX, ThT, pde, n, uhGlobal0{n-1});
		end
		uhSmoGlobal = fem.applyDirichletCondition(ThX, ThT, uhSmo, pde, n);
		for newtonIter = 1: nMaxIt
			NUhSmo = fem.computeNonlinearRHS(ThX, ThT, pde, n, uhSmoGlobal);
			uhS = K\(b + NUhSmo);
			if(norm(uhS - uhSmo) < tol)
				uh = uhS;
				break
			else
				uhSmo = uhS;
				uhSmoGlobal = fem.applyDirichletCondition(ThX, ThT, uhS, pde, n);
			end
		end
		if(newtonIter == nMaxIt)
			warning('The nonlinear solver did not converge');
		end
		uhGlobal0{n} = fem.applyDirichletCondition(ThX, ThT, uh, pde, n);
	end
	% -------------- LOOP OVER DELTA VALUES -------------------- %
	for s = 1 : nMeshRefinements
		pde.delta_ = deltaArray(s);
		for n = 1 : nSteps
			if(n == 1)
				K = fem.computeStiffnessMatrix(ThX, ThT, pde, n);
				b = fem.computeRHS(ThX, ThT, pde, n);
				uhSmo = K\b;
			else
				uhSmo = uh;
				if(n == 2)
					K = fem.computeStiffnessMatrix(ThX, ThT, pde, n);
					K = decomposition(K);
				end
				b = fem.computeRHS(ThX, ThT, pde, n, uhGlobal);
			end
			uhSmoGlobal = fem.applyDirichletCondition(ThX, ThT, uhSmo, pde, n);
			for newtonIter = 1: nMaxIt
				NUhSmo = fem.computeNonlinearRHS(ThX, ThT, pde, n, uhSmoGlobal);
				uhS = K\(b + NUhSmo);
				if(norm(uhS - uhSmo) < tol)
					uh = uhS;
					break
				else
					uhSmo = uhS;
					uhSmoGlobal = fem.applyDirichletCondition(ThX, ThT, uhS, pde, n);
				end
			end
			if(newtonIter == nMaxIt)
				warning('The nonlinear solver did not converge');
			end
			uhGlobal = fem.applyDirichletCondition(ThX, ThT, uh, pde, n);
			aux = fem.getErrorDT_LinfDelta(ThX, ThT, uhGlobal, uhGlobal0{n}, n);
			errDT_Linf(s, p) = max(errDT_Linf(s, p), aux);
			aux = fem.getErrorGradLinfDelta(ThX, ThT, uhGlobal, uhGlobal0{n}, n);
			errGrad_Linf(s, p) = max(errGrad_Linf(s, p), aux);
		end % end of loop over time steps
	end % end of loop over meshes
end % end of loop over degrees

%%---------------- PLOT RATES OF CONVERGENCE --------------------%

FolderName = sprintf('Experiments/delta-convergence/%s/',testName);
mkdir(FolderName);

figure
for p = 1 : pMax
	plotErrors(deltaArray, errDT_Linf(:, p), p);
	hold on;
end
ylabel('$\|\partial_t (u_{h, \tau}^{(\delta)} - u_{h,\tau}^{(0)})\|_{L^{\infty}(0, T; L^2(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\delta$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(deltaArray))
saveas(gcf, strcat(FolderName,'/steepeningDtLinfDelta.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/steepeningDtLinfDelta.eps'), 'epsc')

figure
for p = 1 : pMax
	plotErrors(deltaArray, errGrad_Linf(:, p), p);
	hold on;
end
ylabel('$\|\nabla (u_{h, \tau}^{(\delta)} - u_{h,\tau}^{(0)})\|_{L^{\infty}(0, T; L^2(\Omega)^2)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\delta$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(deltaArray))
saveas(gcf, strcat(FolderName,'/steepeningGradLinfDelta.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/steepeningGradLinfDelta.eps'), 'epsc')
% -----------------------------------------------------------------------------
%																	END OF FILE
% -----------------------------------------------------------------------------