function [hmax] = getHmax(G)
% -----------------------------------------------------------------------------
%  hmax = getHmax() returns the largest cell diameter. 
% -----------------------------------------------------------------------------
	hmax = 0.0;
	numCells = size(G.x_,1)-1;
	for k = 1:numCells,
			dk = G.x_(k+1) - G.x_(k);
			if dk > hmax, hmax = dk; end;
	end;
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
