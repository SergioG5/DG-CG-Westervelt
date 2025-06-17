function [t] = inverseMapPoints(T,x)
% ------------------------------------------------------------------------------
% Get a list of mapped points of the refence 1D-simplex. 
%
% INPUT:
%        t is a N x 1 array of coordinates of points in an arbitrary simplex.
% OUTPUT:
%      t(k) is the mapping of point x(k) with inverse local mapping T<--- To.
% ------------------------------------------------------------------------------
	if size(x,2) ~= 1
		error('Hmmm! in Simplex1D::mapPoints() input must be a Nx1 array');
	end
	t = (2*x-T.v_(2) - T.v_(1))/(T.v_(2)-T.v_(1));
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