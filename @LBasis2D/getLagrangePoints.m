function lagrangePoints = getLagrangePoints(B)
	p = B.deg_;
	N = (p+1)*(p+2)/2;
	lagrangePoints = zeros(N, 2);
	l1 = 1;
	for i = 1: p + 1
		yi = (i-1)/p;
		npi = p + 2 - i;
		xValues = linspace(0, 1.0 - yi, npi);
		yValues = yi*ones(npi, 1);
		lagrangePoints(l1: l1 + npi - 1, 1) = xValues;
		lagrangePoints(l1: l1 + npi - 1, 2) = yValues;
		l1 = l1 + npi;
	end
	lagrangePoints = lagrangePoints(B.map_, :);
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