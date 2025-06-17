function  [] = getSimplex(G, kCell, T)
% -----------------------------------------------------------------------------
% Sets the geometry of cell kCell in Simplex3D T
% -----------------------------------------------------------------------------
	T.setVertices(G.nodes_.xList_(:,G.cells_.vList_(:,kCell)), ...
		G.cells_.vList_(:, kCell), G.cells_.eList_(:, kCell));
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