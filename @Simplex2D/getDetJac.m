function [detJ] = getDetJac(T)
% ------------------------------------------------------------------------------
% area = getArea() gets the area of the simplex.
% ------------------------------------------------------------------------------
	detJ = det(T.j_);
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