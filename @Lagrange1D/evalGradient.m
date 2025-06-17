function [dLdx] = evalGradient(B,X)
% -----------------------------------------------------------------------------
% [dLdx] = evalGradient(X) evaluates the gradient (derivative) of all the
% Lagrange basis functions at each point in the array X. 
%
%  Input:
%
%     X is a N x 1 array of points in [-1 1]
%
% Output:
%
%   dLdx is a dim x N matrix, such that dLdx(i,j) = d(Li)/dx(Xj),
%
%        where Li is the i-th function of the basis.
% -----------------------------------------------------------------------------

	[N, dim2] = size(X);
	if dim2 ~= 1
		error('Dimensions of input array should be N x 1');
	end
	dLegdx = B.Leg_.evalGradient(X);
	dLdx = B.LLM_*dLegdx;
	% convert to zero small entries
	Sp = abs(dLdx) < 5*eps; dLdx(Sp) = 0.0;

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