function [node] = getNode(G, kNode)
% ------------------------------------------------------------------------------
%  G.getNode(k) returns a 2 x 1 array with the coordinates of node k.
% ------------------------------------------------------------------------------
	node = G.nodes_.xList_(:,kNode);
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