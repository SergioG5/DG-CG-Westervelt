function [] = getCell(G, kCell, cell)
% -----------------------------------------------------------------------------
% Returns connectivity of cell k on the object cell of class DGCell2D
%
% -----------------------------------------------------------------------------
	cell.eList_ = G.dg_.eList_(:,kCell);
	cell.nList_ = G.dg_.nList_(:,kCell);
	cell.sList_ = G.dg_.sList_(:,kCell);
	cell.oList_ = G.dg_.rList_(:,kCell);
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