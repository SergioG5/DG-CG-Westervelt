classdef Legendre1D
% -----------------------------------------------------------------------------
% Legendre1D   Legendre orthogonal basis, in the reference 1D-simplex [0 1].
%              Maximum degree allowed is 10.
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Access = private )
	deg_;   % maximum degree
	dim_;   % dimension of polynomial space
end

% -----------------------------------------------------------------------------
%                            METHODS SECTION
% -----------------------------------------------------------------------------
methods

function B = Legendre1D(p)
% -----------------------------------------------------------------------------
% Constructor  
% 
%     B = Legendre1D(p) where polynomial degree 0 <= p <= 10.
% -----------------------------------------------------------------------------

	if int8(p) < 0
		error('Hmmm! in Legendre1D() input parameter should be an integer > 0');
	end
	B.deg_ = p;
	B.dim_ = p + 1;
end

	[dim] = getDimension(B);
	[deg] = getDegree(B);
	LX = evalBasis(B, x);

end % ------------------ END OF METHOD SECTION ------------------------------

% -----------------------------------------------------------------------------
%                          END OF CLASS DEFINITION
% -----------------------------------------------------------------------------
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