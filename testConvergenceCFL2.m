clear all
close all
clc
% -------------------------------------------------------------
% 					TEST PARAMETERS 					
% -------------------------------------------------------------
testNumber = 1;
p = 1;
qMin = 2;
qMax = 4;
nTimeRefinements = 4;
ICopt = 'Ritz';
BCopt = 'Projection';
format short e;

tol = 1.0e-12;
nMaxIt = 15;

% --------------------- MESHES ----------------------------- %

GRID = {'grid-5x5.mat';'grid-10x10.mat';'grid-20x20.mat'; 
'grid-40x40.mat';'grid-80x80.mat'};
N = [50, 100, 200, 400];

% ------------------------------------------------------------- %
% 					ERROR VECTORS  					
% ------------------------------------------------------------- %
nDOFS = zeros(nTimeRefinements, qMax);
errL2 = zeros(nTimeRefinements, qMax);
errDT = zeros(nTimeRefinements, qMax);
errDT_L2T = zeros(nTimeRefinements, qMax);
errDT_Linf = zeros(nTimeRefinements, qMax);
errGrad_Linf = zeros(nTimeRefinements, qMax);
errL2T = zeros(nTimeRefinements, qMax);
errL2H1 = zeros(nTimeRefinements, qMax);
ht = zeros(nTimeRefinements, 1);

% ------------------------------------------------------------- %
% 					TESTS DESCRIPTION  					
% ------------------------------------------------------------- % 
switch(testNumber)
	case 1
		testName = "Gomez-Meliani";
        pde.c_ = 1.0;
		pde.delta_ = 6.0e-9;
		pde.k_ = 0.5;
		A = 1.0e-2;
		omega = 4.5*pi;
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
				A*pde.k_*omega^2*sin(pi*x).*sin(pi*y) ...
				+ 2*A*pde.k_*omega^2*sin(pi*x).*sin(pi*y).*cos(omega*t).^2);
		pde.tf_ = 1.0;
end

% ------------------------------------------------------------- %
% 					MAIN LOOP  					
% ------------------------------------------------------------- %

% ------------------LOOP OVER DEGREES --------------------- %
for q = qMin : qMax
    fem = FEM2D(p, q, ICopt, BCopt);
	
	% -------------- LOOP OVER MESHES -------------------- %
	for s = 1 : nTimeRefinements
		Grid= strcat('Grids/Dirichlet/',GRID{1});
        ThX = Grid2D(Grid);
		tPoints = linspace(0, pde.tf_, N(s) + 1);
		ht(s) = pde.tf_/N(s);
		ThT = Grid1D(tPoints');
		nSteps = N(s);
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
			aux = fem.getErrorL2(ThX, ThT, uhGlobal, pde.u_, n);
			errL2(s, q) = errL2(s, q) + aux^2;
			aux = fem.getErrorGrad(ThX, ThT, uhGlobal, pde.ux_, pde.uy_, n);
			errL2H1(s, q) = errL2H1(s, q) + aux^2;
			aux = fem.getErrorDT(ThX, ThT, uhGlobal, pde.ut_, n);
			errDT(s, q) = errDT(s, q) + aux^2;
			aux = fem.getErrorDT_Linf(ThX, ThT, uhGlobal, pde.ut_, n);
			errDT_Linf(s, q) = max(errDT_Linf(s, q), aux);
			aux = fem.getErrorGradLinf(ThX, ThT, uhGlobal, pde.ux_, pde.uy_, n);
			errGrad_Linf(s, q) = max(errGrad_Linf(s, q), aux);
		end % end of loop over time steps
		errL2T(s, q) = fem.getErrorL2T(ThX, ThT, uhGlobal, pde.u_, n);
		errDT_L2T(s, q) = fem.getErrorDT_L2T(ThX, ThT, uhGlobal, pde.ut_, n);
	end % end of loop over meshes
end % end of loop over degrees
errL2 = sqrt(errL2);
errL2H1 = sqrt(errL2H1);
errDT = sqrt(errDT);

%%---------------- PLOT RATES OF CONVERGENCE --------------------%

FolderName = sprintf('Experiments/CFL2/%s/',testName);
mkdir(FolderName);

figure
for q = qMin : qMax
	plotErrorsTau(ht, errL2(:, q), q);
	hold on;
end
ylabel('$\|u - u_{h, \tau}\|_{L^2(Q_T)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothL2tau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2tau.eps'), 'epsc')


figure
for q = qMin : qMax
	plotErrorsTau(ht, errL2H1(:, q), q);
	hold on;
end
ylabel('$\|u - u_{h, \tau}\|_{L^2(0, T; H^1(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothL2H1tau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2H1tau.eps'), 'epsc')

figure
for q = qMin : qMax
	plotErrorsTau(ht, errL2T(:, q), q);
	hold on;
end
ylabel('$\|u - u_{h, \tau}\|_{L^2(\mathcal{F}^T)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothL2Ftau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2Ftau.eps'), 'epsc')

figure
for q = qMin : qMax
	plotErrorsTau(ht, errDT_L2T(:, q), q);
	hold on;
end
ylabel('$\|\partial_t (u - u_{h, \tau})\|_{L^2(\mathcal{F}^T))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothDtFtau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtFtau.eps'), 'epsc')

figure
for q = qMin : qMax
	plotErrorsTau(ht, errDT(:, q), q);
	hold on;
end
ylabel('$\|\partial_t (u - u_{h, \tau})\|_{L^2(Q_T))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothDtL2tau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtL2tau.eps'), 'epsc')

figure
for q = qMin : qMax
	plotErrorsTau(ht, errDT_Linf(:, q), q);
	hold on;
end
ylabel('$\|\partial_t (u - u_{h, \tau})\|_{L^{\infty}(0, T; L^2(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothDtLinftau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtLinftau.eps'), 'epsc')

figure
for q = qMin : qMax
	plotErrorsTau(ht, errGrad_Linf(:, q), q);
	hold on;
end
ylabel('$\|\nabla (u - u_{h, \tau})\|_{L^{\infty}(0, T; L^2(\Omega)^2)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$\tau$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(ht))
xlim([ht(nTimeRefinements) - 5.0e-3, ht(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothGradLinftau.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothGradLinftau.eps'), 'epsc')
% -----------------------------------------------------------------------------
%																	END OF FILE
% -----------------------------------------------------------------------------