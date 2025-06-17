function [jac] = computeJacobians(G)
% -----------------------------------------------------------------------------
%  jac = getJacobians() returns the list of cell's jacobian relative to 
%                       reference Simplex1D. 
% -----------------------------------------------------------------------------
	numCells = size(G.x_,1)-1;
	jac = zeros(numCells,1);
	T = Simplex1D();
	for k = 1:numCells
		T.setVertices(G.x_(k),G.x_(k+1));
		jac(k) = T.getJacobian();
	end
end



% ------------------------------------------------------------------------------
%                               END OF FILE
% ------------------------------------------------------------------------------




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