classdef Grid1D < handle

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Access = private )
	x_;     % list of nodes
	jac_;	% list of jacobians
end

% -----------------------------------------------------------------------------
%                          METHODS SECTION
% -----------------------------------------------------------------------------
methods

function G = Grid1D(x)
% ----------------------------------------------------------------------------
% Constructor:
%
%   G = Grid1D(x) creates a one dimensional grid with node coordinates given 
%                 in the input Nx1 array x.
% ----------------------------------------------------------------------------
	if size(x,2) ~= 1
		error('Hmmm! input array must be of size Nx1'); 
	end
	G.x_ = x;
	G.jac_ = G.computeJacobians;
end

			[x] 		= getNodes(G);
			[numCells] 	= getNumCells(G);
			[] 			= getGeometry(G,k,T);
			[hmin]		= getHmin(G);
			[hmax] 		= getHmax(G);
			[jac] 		= getJacobians(G);
			[] 			= plot(G);

end												% END OF METHOD SECTION

% -----------------------------------------------------------------------------
%                          END OF CLASS DEFINITION
% -----------------------------------------------------------------------------
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