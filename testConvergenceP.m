% profile on
clear all
close all
clc
% -------------------------------------------------------------
% 					TEST PARAMETERS 					
% -------------------------------------------------------------
testNumber = 1;
pMin       = 2;
pMax       = 9;
nMeshRefinements = 1;
ICopt = 'Ritz';
BCopt = 'Projection';
format short e;

tol = 1.0e-10;
nMaxIt = 1000;

% --------------------- MESHES ----------------------------- %

GRID = {'grid-5x5.mat';'grid-10x10.mat';'grid-20x20.mat'; 
'grid-40x40.mat';'grid-80x80.mat'};
N = [5, 10, 20, 40, 80];

% ------------------------------------------------------------- %
% 					ERROR VECTORS  					
% ------------------------------------------------------------- %
nDOFS = zeros(pMax, 1);
errL2 = zeros(pMax, 1);
errDT = zeros(pMax, 1);
errDT_L2T = zeros(pMax, 1);
errDT_Linf = zeros(pMax, 1);
errGrad_Linf = zeros(pMax, 1);
errL2T = zeros(pMax, 1);
errL2H1 = zeros(pMax, 1);

% ------------------------------------------------------------- %
% 					TESTS DESCRIPTION  					
% ------------------------------------------------------------- % 
switch(testNumber)
	case 1
		testName = "Gomez-Meliani-oscillatory";
        pde.c_ = 1.0;
		pde.delta_ = 6.0e-9;
		pde.k_ = 0.5;
		A = 1.0e-2;
		omega = 9*pi/2;
		pde.u_ = @(x, y, t) A*sin(omega*t).*sin(pi*x).*sin(pi*y);
        pde.ut_ = @(x, y, t) A*omega*cos(omega*t).*sin(pi*x).*sin(pi*y);
		pde.ux_ = @(x, y, t) A*pi*sin(omega*t).*cos(pi*x).*sin(pi*y);
		pde.uy_ = @(x, y, t) A*pi*sin(omega*t).*cos(pi*y).*sin(pi*x);
		pde.u0_ = @(x, y) pde.u_(x, y, 0.0);
		pde.u0x_ = @(x, y) pde.ux_(x, y, 0.0);
		pde.u0y_ = @(x, y) pde.uy_(x, y, 0.0);
        pde.v0_ = @(x, y) pde.ut_(x, y, 0.0);
		pde.f_ = @(x, y, t) A*sin(pi*x).*sin(pi*y).*(2*pde.c_^2*pi^2*sin(omega*t) ...
				- omega^2*sin(omega*t) + 2*pde.delta_*omega*pi^2*cos(omega*t) - ...
				A*pde.k_*omega^2*sin(pi*x).*sin(pi*y) + 2*A*pde.k_*omega^2*sin(pi*x).*sin(pi*y).*cos(omega*t).^2);
		pde.tf_ = 1.0;
end

% ------------------------------------------------------------- %
% 					MAIN LOOP  					
% ------------------------------------------------------------- %

% ------------------LOOP OVER DEGREES --------------------- %
for p = pMin : pMax
	dim = (p + 1)*(p + 2)/2;
    fem = FEM2D(p, p, ICopt, BCopt);
	
	% -------------- LOOP OVER MESHES -------------------- %
		Grid= strcat('Grids/Dirichlet/',GRID{1});
        ThX = Grid2D(Grid);
		tPoints = linspace(0, pde.tf_, N(1) + 1);
		ht = pde.tf_/N(1);
		ThT = Grid1D(tPoints');
		nSteps = N(1);
		for n = 1 : nSteps
			if(n == 1)
				K = fem.computeStiffnessMatrix(ThX, ThT, pde, n);
				b = fem.computeRHS(ThX, ThT, pde, n);
				uhSmo = K\b;
				nDOFS(p) = length(uhSmo)*nSteps;
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
			aux = fem.getErrorL2(ThX, ThT, uhGlobal, pde.u_, n);
			errL2(p) = errL2(p) + aux^2;
			aux = fem.getErrorGrad(ThX, ThT, uhGlobal, pde.ux_, pde.uy_, n);
			errL2H1(p) = errL2H1(p) + aux^2;
			aux = fem.getErrorDT(ThX, ThT, uhGlobal, pde.ut_, n);
			errDT(p) = errDT(p) + aux^2;
			aux = fem.getErrorDT_Linf(ThX, ThT, uhGlobal, pde.ut_, n);
			errDT_Linf(p) = max(errDT_Linf(p), aux);
			aux = fem.getErrorGradLinf(ThX, ThT, uhGlobal, pde.ux_, pde.uy_, n);
			errGrad_Linf(p) = max(errGrad_Linf(p), aux);
		end % end of loop over time steps
		errL2T(p) = fem.getErrorL2T(ThX, ThT, uhGlobal, pde.u_, n);
		errDT_L2T(p) = fem.getErrorDT_L2T(ThX, ThT, uhGlobal, pde.ut_, n);
end % end of loop over degrees
errL2 = sqrt(errL2);
errL2H1 = sqrt(errL2H1);
errDT = sqrt(errDT);

%%---------------- PLOT RATES OF CONVERGENCE --------------------%
FolderName = sprintf('Experiments/p-convergence/%s/',testName);
mkdir(FolderName);

figure
plotErrorsP(nthroot(nDOFS, 3), errL2, 1);
ylabel('$\|u - u_{h,\tau}\|_{L^2(Q_T)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothL2P.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2P.eps'), 'epsc')


figure
plotErrorsP(nthroot(nDOFS, 3), errL2H1, 1);
ylabel('$\|u - u_{h,\tau}\|_{L^2(0, T; H^1(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothL2H1P.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2H1P.eps'), 'epsc')

figure
plotErrorsP(nthroot(nDOFS, 3), errL2T, 1);
ylabel('$\|u - u_{h,\tau}\|_{L^2(\mathcal{F}^T)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothL2FP.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2FP.eps'), 'epsc')

figure
plotErrorsP(nthroot(nDOFS, 3), errDT_L2T, 1);
ylabel('$\|\partial_t (u - u_{h,\tau})\|_{L^2(\mathcal{F}^T))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothDtL2FP.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtL2FP.eps'), 'epsc')

figure
plotErrorsP(nthroot(nDOFS, 3), errDT, 1);
ylabel('$\|\partial_t (u - u_{h,\tau})\|_{L^2(Q_T))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothDtL2P.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtL2P.eps'), 'epsc')

figure
plotErrorsP(nthroot(nDOFS, 3), errDT_Linf, 1);
ylabel('$\|\partial_t (u - u_{h,\tau})\|_{L^{\infty}(0, T; L^2(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothDtP.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtP.eps'), 'epsc')

figure
plotErrorsP(nthroot(nDOFS, 3), errGrad_Linf, 1);
ylabel('$\|\nabla (u - u_{h,\tau})\|_{L^{\infty}(0, T; L^2(\Omega)^2)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\sqrt[3]{N_{\mathrm{DoFs}}}$', 'Interpreter', 'LaTex', 'FontSize', 18);
saveas(gcf, strcat(FolderName,'/smoothGradP.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothGradP.eps'), 'epsc')

% -----------------------------------------------------------------------------
%																	END OF FILE
% -----------------------------------------------------------------------------