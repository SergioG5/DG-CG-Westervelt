function [hmin] = getHmin(G)
% -----------------------------------------------------------------------------
%  hmax = getHmin() returns the smallest cell diameter. 
% -----------------------------------------------------------------------------
	hmin = inf('double');
	numCells = size(G.x_,1)-1;
	for k = 1:numCells
			dk = G.x_(k+1) - G.x_(k);
			if dk < hmin, hmin = dk; end
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
