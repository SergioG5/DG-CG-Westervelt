function [x] = mapPoints(T,xo)
% ------------------------------------------------------------------------------
% x = mapPoints(xo) gets the coordinates of the transformation of points xo
%                   (N x 2) of To by the local mapping To ---> T.
% ------------------------------------------------------------------------------
	if size(xo,2) ~= 2
		error('Hmmm! in Simplex2D::mapPoints() input must be a Nx2 array');
	end
	x = T.j_*xo';
	x(1,:) = x(1,:) + T.v_(1,1);
	x(2,:) = x(2,:) + T.v_(2,1);
	x = x';
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