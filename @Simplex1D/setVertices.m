function [] = setVertices(T,varargin)
% ------------------------------------------------------------------------------
% Set 1D-simplex vertices.
% INPUT:
%   setVertices(v) v is a 1x2 array with vertex coordinates.
%   setVertices(v1,v2) vi coordinates of vertex i
% ------------------------------------------------------------------------------
	switch size(varargin,2)
	case 1
		v1 = varargin{1}(1); v2 = varargin{1}(2);
	case 2
		v1 = varargin{1}; v2 = varargin{2};
	otherwise
	end
	if v1 > v2
		error('Hmmm! wrong vertex order');
	else
		T.v_(1) = v1; T.v_(2) = v2;
			% --------------------- restore jacobian ------------------------------
		T.d_ = T.v_(2) - T.v_(1);
		T.j_ = T.d_;
	end
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