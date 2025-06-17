function [QM] = evalQualityMeasures(G)
	Tk = Simplex2D();
	numCells = G.cells_.numCells_;
	QM.q = zeros(2,numCells);
	for k = 1:numCells
		G.getSimplex(k,Tk);
		m = Tk.getQualityMeasures();
		QM.q(:,k) = m(3:4);
	end
end

% -----------------------------------------------------------------------------
% Created by 
%
% Paul Castillo, paul.castillo@upr.edu
% Department of Mathematical Sciences 
% University of Puerto Rico, Mayaguez Campus (UPRM)
%
% Sergio Gomez, sergio.gomezmacias@unimib.it
% Department of Mathematics and Applications
% University of Milano-Bicocca (UNIMIB)
%
%                                   (2020)
% -----------------------------------------------------------------------------