function [J] = getJacobian(T)
% ------------------------------------------------------------------------------
% Get jacobian of local mapping: refence 1D-simplex to current simplex.
% OUTPUT:
%        jacobian is a 2x2 matrix.
% ------------------------------------------------------------------------------
	J = T.j_;
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