classdef LBasis2D
% -----------------------------------------------------------------------------
% Class definition of a Lagrange interpolatory polynomial basis defined on the 
% reference 2D-simplex, To:  v1 = (0.0), v2 = (1,0) v3 = (0,1). The basis uses
% a uniform node distribution and the maximum degree allowed is 10.
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Constant )
	ZERO_TOL = 1.0e-15;
end

properties ( Access = private )
	map_;
	deg_;
	dim_;
	S_;
	dS_;
end

% -----------------------------------------------------------------------------
%                          PUBLIC METHODS SECTION
% -----------------------------------------------------------------------------

methods

function B = LBasis2D(p)
% -----------------------------------------------------------------------------
% Constructor:
%
%   [B] = LBasis2D(p) creates an Lagrange interpolatory basis object of maximum
%   degree 1 <= p <= 10. The basis uses a uniform node distribution on the 
%   reference 2D-simplex.
% -----------------------------------------------------------------------------
    if ((int8(p) < 1) || (int8(p) > 10))
		error('Hmmm! in LBasis2D() polynomial degree is out of range: [1 .. 10]');
    end
	B.deg_ = p;
	B.dim_ = (p+1)*(p+2)/2;
	B.S_ = ones(p + 1, 3);
	B.dS_ = zeros(p + 1, 3);
	B.map_ = B.getMap;
end

	lagrangePoints = getLagrangePoints(B);
	map = getMap(B);
	[deg] = getDegree(B);
	[dim] = getDimension(B);
	[Px] = evalBasis(B, Xref);
	[dPdx,dPdy] = evalGradient(B,Xref);
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