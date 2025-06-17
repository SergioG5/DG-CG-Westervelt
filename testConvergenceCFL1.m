clear all
close all
clc
% -------------------------------------------------------------
% 					TEST PARAMETERS 					
% -------------------------------------------------------------
testNumber = 1;
pMin       = 1;
pMax       = 3;
q = 2;
nMeshRefinements = 4;
ICopt = 'Ritz';
BCopt = 'Projection';
format short e;

tol = 1.0e-12;
nMaxIt = 15;

% --------------------- MESHES ----------------------------- %

GRID = {'grid-5x5.mat';'grid-10x10.mat';'grid-20x20.mat'; 
'grid-40x40.mat';'grid-80x80.mat'};
N = [5, 5, 5, 5, 5];

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
h = zeros(nMeshRefinements, 1);
hx = zeros(nMeshRefinements, 1);

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
		omega = pi/3;
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
		pde.tf_ = 100.0;
end

% ------------------------------------------------------------- %
% 					MAIN LOOP  					
% ------------------------------------------------------------- %

% ------------------LOOP OVER DEGREES --------------------- %
for p = pMin : pMax
	dim = (p + 1)*(p + 2)/2;
    fem = FEM2D(p, q, ICopt, BCopt);
	
	% -------------- LOOP OVER MESHES -------------------- %
	for s = 1 : nMeshRefinements
		Grid= strcat('Grids/Dirichlet/',GRID{s});
        ThX = Grid2D(Grid);
		tPoints = linspace(0, pde.tf_, N(s) + 1);
		ht = pde.tf_/N(s);
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
			[~ ,hx(s)] = ThX.getGridSize();
			h(s) = hx(s);
			aux = fem.getErrorL2(ThX, ThT, uhGlobal, pde.u_, n);
			errL2(s, p) = errL2(s, p) + aux^2;
			aux = fem.getErrorGrad(ThX, ThT, uhGlobal, pde.ux_, pde.uy_, n);
			errL2H1(s, p) = errL2H1(s, p) + aux^2;
			aux = fem.getErrorDT(ThX, ThT, uhGlobal, pde.ut_, n);
			errDT(s, p) = errDT(s, p) + aux^2;
			aux = fem.getErrorDT_Linf(ThX, ThT, uhGlobal, pde.ut_, n);
			errDT_Linf(s, p) = max(errDT_Linf(s, p), aux);
			aux = fem.getErrorGradLinf(ThX, ThT, uhGlobal, pde.ux_, pde.uy_, n);
			errGrad_Linf(s, p) = max(errGrad_Linf(s, p), aux);
		end % end of loop over time steps
		errL2T(s, p) = fem.getErrorL2T(ThX, ThT, uhGlobal, pde.u_, n);
		errDT_L2T(s, p) = fem.getErrorDT_L2T(ThX, ThT, uhGlobal, pde.ut_, n);
	end % end of loop over meshes
end % end of loop over degrees
errL2 = sqrt(errL2);
errL2H1 = sqrt(errL2H1);
errDT = sqrt(errDT);

%%---------------- PLOT RATES OF CONVERGENCE --------------------%

FolderName = sprintf('Experiments/CFL/%s/',testName);
mkdir(FolderName);

figure
for p = 1 : pMax
	plotErrors(h, errL2(:, p), p);
	hold on;
end
ylabel('$\|u - u_h\|_{L^2(Q_T)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothL2H.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2H.eps'), 'epsc')


figure
for p = 1 : pMax
	plotErrors(h, errL2H1(:, p), p);
	hold on;
end
ylabel('$\|u - u_h\|_{L^2(0, T; H^1(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothL2H1H.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2H1H.eps'), 'epsc')

figure
for p = 1 : pMax
	plotErrors(h, errL2T(:, p), p);
	hold on;
end
ylabel('$\|u - u_h\|_{L^2(\mathcal{F}^T)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothL2FH.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothL2FH.eps'), 'epsc')

figure
for p = 1 : pMax
	plotErrors(h, errDT_L2T(:, p), p);
	hold on;
end
ylabel('$\|\partial_t (u - u_h)\|_{L^2(\mathcal{F}^T))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothDtFH.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtFH.eps'), 'epsc')

figure
for p = 1 : pMax
	plotErrors(h, errDT(:, p), p);
	hold on;
end
ylabel('$\|\partial_t (u - u_h)\|_{L^2(Q_T))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothDtL2H.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtL2H.eps'), 'epsc')

figure
for p = 1 : pMax
	plotErrors(h, errDT_Linf(:, p), p);
	hold on;
end
ylabel('$\|\partial_t (u - u_{h,\tau})\|_{L^{\infty}(0, T; L^2(\Omega))}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothDtLinfH.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothDtLinfH.eps'), 'epsc')

figure
for p = 1 : pMax
	plotErrors(h, errGrad_Linf(:, p), p);
	hold on;
end
ylabel('$\|\nabla (u - u_{h, \tau})\|_{L^{\infty}(0, T; L^2(\Omega)^2)}$',...%
			'Interpreter','LaTex', 'FontSize', 18);
xlabel('$h$', 'Interpreter', 'LaTex', 'FontSize', 18);
xtickformat('%.1e')
xticks(flip(h))
xlim([h(nMeshRefinements) - 5.0e-3, h(1) + 1.0e-2]);
saveas(gcf, strcat(FolderName,'/smoothGradLinfH.fig'), 'fig')
saveas(gcf, strcat(FolderName,'/smoothGradLinfH.eps'), 'epsc')
% ------------------------------------------------------------------------
%															END OF FILE
% ------------------------------------------------------------------------