function [numEdges] = getNumEdges(G)
% -----------------------------------------------------------------------------
% Returns the total number of physical edges: interior and boundary.
% -----------------------------------------------------------------------------
	numEdges = G.edges_.numEdges_;
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