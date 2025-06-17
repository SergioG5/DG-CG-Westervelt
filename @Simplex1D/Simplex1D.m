classdef Simplex1D < handle
% -----------------------------------------------------------------------------
% 1D-Simplex class definition. The reference 1D-simplex is the interval [0 1]
% -----------------------------------------------------------------------------

% -----------------------------------------------------------------------------
%                          PROPERTIES SECTION
% -----------------------------------------------------------------------------

properties ( Access = private )
	v_; % 1x2 array of vertex coordinates
	j_; % Jacobian of linear mapping
	d_; % simplex diameter
end

% -----------------------------------------------------------------------------
%                          METHODS SECTION
% -----------------------------------------------------------------------------
methods

function T = Simplex1D(varargin)
% -----------------------------------------------------------------------------
% 1D-simplex constructor:
%
% INPUT:
%
%     Simplex1D()      the reference 1D-simplex [-1 1] is created.
%     Simplex1D(v)     v is a 1x2 array with the vertices
%     Simplex1D(v1,v2) v1 left end point v2 right end point 
% -----------------------------------------------------------------------------
	T.v_ = zeros(1,2);
	optargin = size(varargin,2);
	switch optargin
	case 0
		T.v_ = [0 1];     % REFERENCE 1D-SIMPLEX
	case 1
		T.v_ = varargin{1};
	case 2
		T.v_(1,1) = varargin{1};
		T.v_(1,2) = varargin{2};
	otherwise
		error('Hmmm! in Simplex1D(): number of input parameters must be 0,1 or 2')
	end
	% ------------------------------- jacobian ----------------------------------
	T.d_ = T.v_(2)-T.v_(1);
	T.j_ = T.d_;
	%------------------------ normals and edge lengths --------------------------

end


			[] 		= setVertices(T,varargin);
	[vertexList] 	= getVertices(T);
			[J]		= getJacobian(T);
			[x] 	= mapPoints(T,xo);
			[t] 	= inverseMapPoints(T,x);
	[barycenter] 	= getBarycenter(T);
		[diam] 		= getDiameter(T);

end											% END OF METHODS SECTION

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