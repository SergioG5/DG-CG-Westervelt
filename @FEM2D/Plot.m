function [gcfLdg] = Plot(fem, ThX, uh)
	n = 100;
	[v, w] = meshgrid(0 : 1/n : 1, 0 : 1/n : 1);
	X = v.*(1 - w);
	Y = w;
	aux = cell(length(Y), 1);
	Z = zeros(size(X));
	Xk = zeros(size(X));
	Yk = zeros(size(X));
	PTTop = fem.op2D_.PTTopQ_;
	PXT = cell(n + 1, 1);
	for i = 1 : n + 1
		PX = fem.op2D_.BX_.evalBasis([X(:,i), Y(:,i)]);
		PXT{i} = kron(PTTop, PX);
	end
	numCells = ThX.getNumCells();
	colormap jet
	shading interp;
	hold on;
	Tk = Simplex2D();
	for kCell = 1:numCells
		ThX.getSimplex(kCell, Tk);
		DOF = fem.getDOF(ThX, kCell);
		for i = 1 : n+1
			aux{i} = Tk.mapPoints([X(:,i), Y(:,i)]);
			Xk(:,i) = aux{i}(:, 1);
			Yk(:,i) = aux{i}(:, 2);
			Z(:,i) = PXT{i}'*uh(DOF.ptrsGlobalDofs);
		end
		gcfLdg = surf(Xk, Yk, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
	end
	grid on;
	box on;
	hold off;
	set(gca, 'LineWidth', 1.0, 'Fontsize', 12.0, 'FontWeight', 'bold');
end
% ------------------------------------------------------------------------------
%                               END OF FILE
% ------------------------------------------------------------------------------


% -----------------------------------------------------------------------------
% Created by 
%
% Sergio Gomez, sergio.gomezmacias@unimib.it
% Department of Mathematics and Applications
% University of Milano-Bicocca (UNIMIB)
%
%                                   (2025)
% -----------------------------------------------------------------------------