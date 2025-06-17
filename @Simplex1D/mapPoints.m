function [x] = mapPoints(T,xo)
% ------------------------------------------------------------------------------
% Get a list of mapped points of the refence 1D-simplex. 
%
% INPUT:
%        xo is a N x 1 array of coordinates of points in the reference simplex.
% OUTPUT:
%      x(k) is the mapping of point xo(k) with local mapping To ---> T.
% ------------------------------------------------------------------------------
	if size(xo,2) ~= 1
		error('Hmmm! in Simplex1D::mapPoints() input must be a Nx1 array');
	end
	x = T.j_*xo + T.v_(1);
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