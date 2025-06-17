function [hmin, hmax] = getGridSize(G)
% -----------------------------------------------------------------------------
% [hmin,hmax] = getGridSize() returns minimum and maximum edge length of G.
% -----------------------------------------------------------------------------
	hmax = 0.0;
	hmin = inf('double');
	numEdges = G.edges_.numEdges_;
	for k = 1:numEdges
			v1 = G.nodes_.xList_(:,G.edges_.vList_(1,k));
			v2 = G.nodes_.xList_(:,G.edges_.vList_(2,k));
			d12 = norm(v1-v2);
			if d12 > hmax, hmax = d12; end
			if d12 < hmin, hmin = d12; end
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