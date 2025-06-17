classdef Lagrange1D
% -----------------------------------------------------------------------------
% Legendre1D   definition of a class for a polynomial basis in the reference 
%              1D-cell [0 1], using Lagrange polynomials.
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Access = private )
	deg_;   % maximum degree
	dim_;   % dimension of polynomial space
	Leg_;   % auxiliary Legendre polynomials
	pts_;   % interpolatory points
	LLM_;   % Legendre to Lagrange matrix.
end


% -----------------------------------------------------------------------------
%                          PUBLIC METHODS SECTION
% -----------------------------------------------------------------------------

methods (Access = public)

function [B] = Lagrange1D(p, varargin)
% -----------------------------------------------------------------------------
% Constructor:
%
%   [B] = Lagrange1D(p) 
%         creates an object for a Lagrange polynomial basis on the reference 
%         cell [0 1]; of maximum degree p and uniformly distributed 
%         interpolation nodes.
%
%   [B] = Lagrange1D(p,nodeList) 
%         creates an object for a Lagrange polynomial basis on the reference 
%         cell [0 1]; of maximum degree p whose interpolation nodes are given 
%         in the (p+1) x 1 input array nodeList.
%
% -----------------------------------------------------------------------------
	if int8(p) < 0
		error('Polynomial degree %d should be greater or equal than 0',p);
	end
	B.deg_ = p;
	B.dim_ = p + 1;
	B.Leg_ = Legendre1D(p);
	switch ( nargin )
	case 1
		B.pts_ = linspace(0.0, 1.0, B.deg_ + 1)';
	case 2
		nodes = sort(varargin{1});
		if ( length(nodes) ~= p+1 )
			error('Incompatible number of interpolation nodes');
		elseif(length(unique(nodes)) ~= p+1)
			error('Repeated interpolation nodes');
		elseif( nodes(1) < 0 ||nodes(p+1) > 1)
				error('Interpolation nodes out of range [-1 1]');
		else
			B.pts_ = nodes;
		end
	otherwise
		error('Invalid number of input parameters in constructor');
	end
	
	PX = B.Leg_.evalBasis(B.pts_);
	B.LLM_ = eye(B.dim_)/PX;

end

	% ------------------------------------------------------------------------
		 [dim] = dimension(B);
		 [deg] = degree(B);
	 [nodes] = interpolationNodes(B);
		  [LX] = evalBasis(B, X);
		[dLdx] = evalGradient(B, X);
        [d2Ldx2] = evalSecondDerivatives(B, X);
	% ------------------------------------------------------------------------

end 									% END OF PUBLIC METHODS SECTION

% -----------------------------------------------------------------------------
%                          END OF CLASS DEFINITION
% -----------------------------------------------------------------------------
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