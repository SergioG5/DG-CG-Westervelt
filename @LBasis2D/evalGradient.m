function [dPdx,dPdy] = evalGradient(B, Xref)
% -----------------------------------------------------------------------------
% [dPdx,dPdy] = evalGradient(X) evaluates the gradient components of all the 
% basis functions for all the points in the [N x 2] array X of points in the 
% reference 2D-simplex. dPdx and dPdy are [dim x N] matrices, where 
%      dPdx(i,j) = d(Pi)/dx(Xj) and dPdy(i,j0 =  d(Pi)/dy(Xj) 
% and Pi is the i-th function of the basis.
% -----------------------------------------------------------------------------

	if size(Xref,2) ~= 2
		error('Hmmm! in LBasis2D::evalGradient() dimension of input array Nx2');
	end
	N = size(Xref,1);
	dPdx = zeros(B.dim_,N);
	dPdy = zeros(B.dim_,N);
	for j = 1 : N
		z2 = Xref(j,1);
		z3 = Xref(j,2);
		z1 = 1.0 - z2 -z3;
		for k = 1 : B.deg_
			a0 = B.deg_/k;
		
			a1 = (B.deg_ * z1 - k + 1.0)/k;
			B.S_(k+1,1)  = a1*B.S_(k,1);
			B.dS_(k+1,1) = a1*B.dS_(k,1) + a0*B.S_(k,1);
		
			a1 = (B.deg_ * z2 - k + 1.0)/k;
			B.S_(k+1,2)  = a1*B.S_(k,2);
			B.dS_(k+1,2) = a1*B.dS_(k,2) + a0*B.S_(k,2);
		
			a1 = (B.deg_ * z3 - k + 1.0)/k;
			B.S_(k+1,3)  = a1*B.S_(k,3);
			B.dS_(k+1,3) = a1*B.dS_(k,3) + a0*B.S_(k,3);
		end

		k = 1;
		for k3 = 0 : B.deg_
			i3 = k3 + 1;
			for k2 = 0 : B.deg_ - k3
				i2 = k2 + 1;
				i1 = B.deg_ - k2 - k3 + 1;
				dPdx(k,j) = B.S_(i3,3)*(B.dS_(i2,2)*B.S_(i1,1) - B.S_(i2,2)*B.dS_(i1,1));
				dPdy(k,j) = B.S_(i2,2)*(B.dS_(i3,3)*B.S_(i1,1) - B.S_(i3,3)*B.dS_(i1,1));
				k = k + 1;
			end
		end
	end
	% convert to zero small entries
	Sp = abs(dPdx) < B.ZERO_TOL; dPdx(Sp) = 0.0;
	Sp = abs(dPdy) < B.ZERO_TOL; dPdy(Sp) = 0.0;
	dPdx = dPdx(B.map_, :);
	dPdy = dPdy(B.map_, :);
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