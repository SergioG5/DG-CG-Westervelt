% profile on
clear all
close all
clc
% -------------------------------------------------------------
% 					TEST PARAMETERS 					
% -------------------------------------------------------------
testNumber = 3;
pMin       = 4;
pMax       = 4;
q = 4;
ICopt = 'Ritz';
BCopt = 'Projection';

tol = 1.0e-12;
nMaxIt = 15;

% --------------------- MESHES ----------------------------- %

GRID = {'grid-10x10.mat'};
N = 20;

% ------------------------------------------------------------- %
% 					TESTS DESCRIPTION  					
% ------------------------------------------------------------- % 
switch(testNumber)
	case 3
		testName = "Steepening-wave-front";
		pde.delta_ = 6.0e-9;
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
uhGlobal0 = cell(N, 1);
for p = pMin : pMax
	dim = (p + 1)*(p + 2)/2;
    fem = FEM2D(p, q, ICopt, BCopt);
	Grid = strcat('Grids/Dirichlet/',GRID{1});
	ThX = Grid2D(Grid);
	tPoints = linspace(0, pde.tf_, N(1) + 1);
	ht = pde.tf_/N(1);
	ThT = Grid1D(tPoints');
	nSteps = N(1);
	% -------------- SOLUTION OF VANISHING LIMIT -------------------- %
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
		figure
		fem.Plot(ThX, uhGlobal0{n});
		view([-45, 25])
		clim([-1.e-8, 2.0e-7])
		zlim([-1.e-8, 2.0e-7])
		colorbar
		xlabel('$x$','Interpreter','LaTex', 'FontSize', 18);
		ylabel('$y$','Interpreter','LaTex', 'FontSize', 18);
		zlabel('$u_{h, \tau}$', 'Interpreter','LaTex', 'FontSize', 18);
	end
end % end of loop over degrees

% -----------------------------------------------------------------------------
%																	END OF FILE
% -----------------------------------------------------------------------------