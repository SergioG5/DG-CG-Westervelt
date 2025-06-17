function [Px] = evalBasis(B,Xref)
% -----------------------------------------------------------------------------
% [Px] = evalBasis(X) evaluates all the basis functions for all the points in 
% the [N x 2] array X of points in the reference 2D-simplex. PX is a [dim x N] 
% matrix, where PX(i,j) = Pi(Xj), Pi is the i-th basis function of the basis.
% -----------------------------------------------------------------------------
	if size(Xref,2) ~= 2
		error('Hmmm! in LBasis2D::evalBasis() dimension of input array Nx2');
	end
	N = size(Xref,1);
	Px = zeros(B.dim_,N);
	for j = 1:N
		z2 = Xref(j,1);
		z3 = Xref(j,2);
		z1 = 1.0 - z2 -z3;
		for k = 1 : B.deg_
			B.S_(k+1,1) = B.S_(k,1)*(B.deg_ * z1 - k + 1.0)/k;
			B.S_(k+1,2) = B.S_(k,2)*(B.deg_ * z2 - k + 1.0)/k;
			B.S_(k+1,3) = B.S_(k,3)*(B.deg_ * z3 - k + 1.0)/k;
		end
		k = 1;
		for k3 = 0 : B.deg_
			for k2 = 0 : B.deg_ - k3
				Px(k,j) = B.S_(B.deg_ - k2 - k3 + 1,1)*B.S_(k2+1,2)*B.S_(k3+1,3);
				k = k + 1;
			end
		end
	end
	% convert to zero small entries
	Sp = abs(Px) < B.ZERO_TOL; Px(Sp) = 0.0;
	Px = Px(B.map_, :);
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