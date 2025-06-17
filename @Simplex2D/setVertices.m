function [] =  setVertices(T,vertexList, vertexIndex, edgeIndex)
% ------------------------------------------------------------------------------
% setVertices(vertexList) sets 2D-simplex vertices whose coordinates are 
%                         given in the input 2x3 array. All internal data, such 
%                         as jacobian, normals, edge lentghs, etc are updated.
% ------------------------------------------------------------------------------
	T.iv_ = vertexIndex;
	T.ie_ = edgeIndex;
	T.v_ = vertexList;
		% --------------------- restore jacobian ------------------------------
	T.j_(:,1) = T.v_(:,2) - T.v_(:,1);
	T.j_(:,2) = T.v_(:,3) - T.v_(:,1);
		% --------------- restore normals and edge lengths --------------------
	% Normal relative to edge 1: (v2 v3)
	nv = [ T.v_(2,3) - T.v_(2,2)
	       T.v_(1,2) - T.v_(1,3) ];
	edgeLength = norm(nv);
	T.e_(1) = edgeLength;
	T.n_(:,1) = nv/edgeLength;
	% Normal relative to edge 2: (v3 v1)
	nv = [ T.v_(2,1) - T.v_(2,3)
	       T.v_(1,3) - T.v_(1,1) ];
	edgeLength = norm(nv);
	T.e_(2) = edgeLength;
	T.n_(:,2) = nv/edgeLength;
	% Normal relative to edge 3: (v1 v2)
	nv = [ T.v_(2,2) - T.v_(2,1)
	       T.v_(1,1) - T.v_(1,2) ];
	edgeLength = norm(nv);
	T.e_(3) = edgeLength;
	T.n_(:,3) = nv/edgeLength;
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