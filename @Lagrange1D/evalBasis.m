function [LX] = evalBasis(B, X)
% -----------------------------------------------------------------------------
% [LX] = evalBasis(X) evaluates all the Lagrange basis functions at each point 
% in the array X. 
%
%  Input:
%
%     X is a N x 1 array of points in [0 1]
%
% Output:
%
%    LX is a dim x N matrix, such that  LX(i,j) = Li(Xj), 
%
%       where Li is the i-th basis function.
% -----------------------------------------------------------------------------

	[N, dim2] = size(X);
	if dim2 ~= 1
		error('Dimensions of input array should be N x 1');
	end
	LegX = B.Leg_.evalBasis(X);
	LX = B.LLM_*LegX;
		% convert to zero small entries
	Sp = abs(LX) < 5*eps; LX(Sp) = 0.0;
	
end % end of function




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