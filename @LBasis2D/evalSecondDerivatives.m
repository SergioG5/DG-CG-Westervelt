function [d2Pdxx, d2Pdyx, d2Pdxy, d2Pdyy] = evalSecondDerivatives(B,Xref)
% -----------------------------------------------------------------------------
% [dPdx,dPdy] = evalGradient(X) evaluates the gradient components of all the 
% basis functions for all the points in the [N x 2] array X of points in the 
% reference 2D-simplex. dPdx and dPdy are [dim x N] matrices, where 
%      dPdx(i,j) = d(Pi)/dx(Xj) and dPdy(i,j0 =  d(Pi)/dy(Xj) 
% and Pi is the i-th function of the basis.
% -----------------------------------------------------------------------------
	p = B.deg_;
	nPoints = size(Xref, 1);
	ONES = ones(1, nPoints);
	ZEROS = zeros(1, nPoints);
	switch(p)
		case 1
			d2Pdxx = zeros(nPoints, 3);
			d2Pdxy = zeros(nPoints, 3);
			d2Pdyx = zeros(nPoints, 3);
			d2Pdyy = zeros(nPoints, 3);
		case 2
			d2Pdxx = zeros(nPoints, 6);
			d2Pdxy = zeros(nPoints, 6);
			d2Pdyx = zeros(nPoints, 6);
			d2Pdyy = zeros(nPoints, 6);
			dPhi1dx = -ONES;
			dPhi1dy = -ONES;
			dPhi2dx = ONES;
			dPhi2dy = ZEROS;
			dPhi3dx = ZEROS;
			dPhi3dy = ONES;
			d2Pdxx(:, 1) = dPhi1dx.*(4.0.*dPhi1dx);
			d2Pdxy(:, 1) = dPhi1dx.*(4.0.*dPhi1dy);
			d2Pdyx(:, 1) = dPhi1dy.*(4.0.*dPhi1dx);
			d2Pdyy(:, 1) = dPhi1dy.*(4.0.*dPhi1dy);
			d2Pdxx(:, 2) = dPhi2dx.*(4.0.*dPhi2dx);
			d2Pdxy(:, 2) = dPhi2dx.*(4.0.*dPhi2dy);
			d2Pdyx(:, 2) = dPhi2dy.*(4.0.*dPhi2dx);
			d2Pdyy(:, 2) = dPhi2dy.*(4.0.*dPhi2dy);
			d2Pdxx(:, 3) = dPhi3dx.*(4.0.*dPhi3dx);
			d2Pdxy(:, 3) = dPhi3dx.*(4.0.*dPhi3dy);
			d2Pdyx(:, 3) = dPhi3dy.*(4.0.*dPhi3dx);
			d2Pdyy(:, 3) = dPhi3dy.*(4.0.*dPhi3dy);
			d2Pdxx(:, 4) = 4.0*(dPhi2dx.*dPhi3dx + dPhi2dx.*dPhi3dx);
			d2Pdxy(:, 4) = 4.0*(dPhi2dx.*dPhi3dy + dPhi2dy.*dPhi3dx);
			d2Pdyx(:, 4) = 4.0*(dPhi2dy.*dPhi3dx + dPhi2dx.*dPhi3dy);
			d2Pdyy(:, 4) = 4.0*(dPhi2dy.*dPhi3dy + dPhi2dy.*dPhi3dy);
			d2Pdxx(:, 5) = 4.0*(dPhi3dx.*dPhi1dx + dPhi3dx.*dPhi1dx);
			d2Pdxy(:, 5) = 4.0*(dPhi3dx.*dPhi1dy + dPhi3dy.*dPhi1dx);
			d2Pdyx(:, 5) = 4.0*(dPhi3dy.*dPhi1dx + dPhi3dx.*dPhi1dy);
			d2Pdyy(:, 5) = 4.0*(dPhi3dy.*dPhi1dy + dPhi3dy.*dPhi1dy);
			d2Pdxx(:, 6) = 4.0*(dPhi1dx.*dPhi2dx + dPhi1dx.*dPhi2dx);
			d2Pdxy(:, 6) = 4.0*(dPhi1dx.*dPhi2dy + dPhi1dy.*dPhi2dx);
			d2Pdyx(:, 6) = 4.0*(dPhi1dy.*dPhi2dx + dPhi1dx.*dPhi2dy);
			d2Pdyy(:, 6) = 4.0*(dPhi1dy.*dPhi2dy + dPhi1dy.*dPhi2dy);
		case 3
			d2Pdxx = zeros(nPoints, 10);
			d2Pdxy = zeros(nPoints, 10);
			d2Pdyx = zeros(nPoints, 10);
			d2Pdyy = zeros(nPoints, 10);
			x = Xref(:, 1);
			y = Xref(:, 2);
			d2Pdxx(:, 1) = 18 - 27*y - 27*x;
			d2Pdxy(:, 1) = 18 - 27*y - 27*x;
			d2Pdyx(:, 1) = 18 - 27*y - 27*x;
			d2Pdyy(:, 1) = 18 - 27*y - 27*x;
			d2Pdxx(:, 2) = 27*x - 9;
			d2Pdxy(:, 2) = ZEROS;
			d2Pdyx(:, 2) = ZEROS;
			d2Pdyy(:, 2) = ZEROS;
			d2Pdxx(:, 3) = ZEROS;
			d2Pdxy(:, 3) = ZEROS;
			d2Pdyx(:, 3) = ZEROS;
			d2Pdyy(:, 3) = 27*y - 9;
			d2Pdxx(:, 4) = -54*y;
			d2Pdxy(:, 4) = 27 - 54*y - 54*x;
			d2Pdyx(:, 4) = 27 - 54*y - 54*x;
			d2Pdyy(:, 4) = -54*x;
			d2Pdxx(:, 5) = 27*y;
			d2Pdxy(:, 5) = 27*x - 9/2;
			d2Pdyx(:, 5) = 27*x - 9/2;
			d2Pdyy(:, 5) = ZEROS;
			d2Pdxx(:, 6) = ZEROS;
			d2Pdxy(:, 6) = 27*y - 9/2;
			d2Pdyx(:, 6) = 27*y - 9/2;
			d2Pdyy(:, 6) = 27*x;
			d2Pdxx(:, 7) = ZEROS;
			d2Pdxy(:, 7) = 9/2 - 27*y;
			d2Pdyx(:, 7) = 9/2 - 27*y;
			d2Pdyy(:, 7) =  36 - 81*y - 27*x;
			d2Pdxx(:, 8) = 27*y;
			d2Pdxy(:, 8) = 27*x + 54*y - 45/2;
			d2Pdyx(:, 8) = 27*x + 54*y - 45/2;
			d2Pdyy(:, 8) = 54*x + 81*y - 45;
			d2Pdxx(:, 9) =  81*x + 54*y - 45;
			d2Pdxy(:, 9) = 54*x + 27*y - 45/2;
			d2Pdyx(:, 9) = 54*x + 27*y - 45/2;
			d2Pdyy(:, 9) = 27*x;
			d2Pdxx(:, 10) =  36 - 27*y - 81*x;
			d2Pdxy(:, 10) = 9/2 - 27*x;
			d2Pdyx(:, 10) = 9/2 - 27*x;
			d2Pdyy(:, 10) = ZEROS;
	end

	d2Pdxx = d2Pdxx';
	d2Pdxy = d2Pdxy';
	d2Pdyx = d2Pdyx';
	d2Pdyy = d2Pdyy';
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