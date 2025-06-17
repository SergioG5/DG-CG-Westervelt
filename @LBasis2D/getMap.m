function map = getMap(B)
	map = zeros(B.dim_, 1);
	p = B.deg_;
	map(1:3) = [1, p + 1, (p+1)*(p+2)/2];
	cont = 4;
	bNode = p + 3;
	for k = 1: p - 2
		dk = p - k - 2;
		map(cont: cont + dk) = bNode: bNode + dk;
		cont = cont + dk + 1;
		bNode = bNode + dk + 3;
	end
	e1 = 2*p + 1;
	for k = 1: p - 1
		dk = p - k;
		map(cont) = e1;
		cont = cont + 1;
		e1 = e1 + dk;
	end
	e2 = p + 2;
	iAux = cont : cont + p - 2;
	for k = 1: p - 1
		dk = p - k + 1;
		map(cont) = e2;
		cont = cont + 1;
		e2 = e2 + dk;
	end
    map(cont:cont + p - 2) = 2: p;
	map(iAux) = flip(map(iAux));
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