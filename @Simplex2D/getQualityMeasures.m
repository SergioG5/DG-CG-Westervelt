function [m] = getQualityMeasures(T)
% ------------------------------------------------------------------------------
% m = getQualityMeasures() gets some quality measures of the simplex.
%
%   m(1): radius of incircle (ir).
%   m(2): radius of circumscribed circle (or).
%   m(3): 2.0*ir/or.
%   m(4): sqrt(3)*ir/diam(T)
% ------------------------------------------------------------------------------
	m = zeros(4,1);
	x = (T.e_(1) + T.e_(2) - T.e_(3))/2.0;
	y = (T.e_(1) + T.e_(3) - T.e_(2))/2.0;
	z = (T.e_(2) + T.e_(3) - T.e_(1))/2.0;
	s  = x+y+z;
	m(1) = sqrt((x*y*z)/s);
	m(2) = (T.e_(1)*T.e_(2)*T.e_(3))/(4.0*s*m(1));
	m(3) = 2.0*m(1)/m(2);
	m(4) = sqrt(3.0)*m(1)/max(T.e_);
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