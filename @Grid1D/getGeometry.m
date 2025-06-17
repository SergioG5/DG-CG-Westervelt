function  [] = getGeometry(G,k,T)
% -----------------------------------------------------------------------------
%  getGeometry(k,T) sets geometry of cell k in Simplex1D T.
% -----------------------------------------------------------------------------
	T.setVertices(G.x_(k), G.x_(k+1));
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