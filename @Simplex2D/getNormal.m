function [ni] = getNormal(T,i)
% ------------------------------------------------------------------------------
% ni = getNormal(T,i) get the outward unitary normal vector at edge i.
% ------------------------------------------------------------------------------
	ni = T.n_(:,i);
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