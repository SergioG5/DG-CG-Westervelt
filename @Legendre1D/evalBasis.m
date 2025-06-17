function LX = evalBasis(B,x)
% -----------------------------------------------------------------------------
% Evaluate basis functions at points in the reference 1D-simplex [0 1].
% 
%  [LX] = evalBasis(X) evaluates all the basis functions at the input nodes X. 
%  X is a [N x 1] array of point coordinates in the reference 1D-simplex [0 1] 
%  LX is a [dim X N] matrix, where LX(i,j) = Li(Xj), Li is the Legendre 
%  polynomial of degree i-1.
% -----------------------------------------------------------------------------
	if size(x,2) ~= 1
		error('Hmmm! in Legendre1D::evalBasis() input should be a N x 1 array');
	end
	LX = ones(length(x),B.deg_+1);
	if B.deg_ > 0
		LX(:, 2) = 2.0*x - 1.0;
		for k = 3 : B.dim_
			a1 = (2*k - 3)/(k - 1);
			a2 = -(k - 2)/(k - 1);
			LX(:,k) = a1 *(2.0*x - 1.0) .* LX(:, k-1) + a2 * LX(:, k-2);
		end
	end
	LX = LX';
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