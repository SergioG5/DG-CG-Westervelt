function [d2Ldx2] = evalSecondDerivatives(B, X)
% -----------------------------------------------------------------------------
% [dLdx] = evalGradient(X) evaluates the gradient (derivative) of all the
% Legendre basis functions for each point in the array X. 
%
%  Input:
%
%     X is a N x 1 array of points in [0 1]
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
	Lx   = ones(N,B.dim_);
	dLdx = zeros(N,B.dim_);
    d2Ldx2 = zeros(N,B.dim_);
	if B.deg_ > 0
		Lx(:, 2) = 2*X - 1;
		dLdx(:,2) = 2.0;
		for k = 3:B.dim_
			a1 = (2*k - 3)/(k - 1);
			a2 = -(k - 2)/(k - 1);
			Lx(:, k) = a1 * (2*X - 1) .* Lx(:,k - 1) + a2 * Lx(:, k - 2);
			dLdx(:, k) = a1 * (2*Lx(:,k-1) + (2*X - 1).* dLdx(:,k-1)) ...
						+ a2 * dLdx(:,k-2);
			d2Ldx2(:, k) = a1*(4*dLdx(:, k - 1) + (2*X - 1).*d2Ldx2(:, k - 1)) ...
						 + a2 * d2Ldx2(:, k - 2);
		end
	end
	d2Ldx2 = d2Ldx2';
	% convert to zero small entries
	Sp = abs(d2Ldx2) < 5*eps; d2Ldx2(Sp) = 0.0;
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