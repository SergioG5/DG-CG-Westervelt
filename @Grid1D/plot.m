function [] = plot(G)
% -----------------------------------------------------------------------------
%  plot() visualizes the current grid. 
% -----------------------------------------------------------------------------
	figure;
	numNodes = size(G.x_,1);
	axis([G.x_(1) G.x_(numNodes) 0 2])
	plot(G.x_,ones(numNodes,1),'o-')
	hold on;
	delta = 0.1;
	% ---------------------------- NODES ---------------------------------------
	y = 1.0 + delta;
	for k = 1:numNodes
		xk = G.x_(k);
		text(xk,y,num2str(k),'Color','b','FontSize',8);
		plot(xk,1,'o','MarkerSize',8,'MarkerEdgeColor','k', ...
				'MarkerFaceColor','b');
	end
	% ---------------------------- CELLS ---------------------------------------
	numCells = numNodes-1;
	gp = (G.x_(2:numNodes) + G.x_(1:numCells))/2.0;
	plot(gp,ones(numCells,1),'o','MarkerSize',6,...
					'MarkerEdgeColor','k','MarkerFaceColor','r');
		y = 1.0 - delta;
	for k = 1:numCells
		text(gp(k),y,num2str(k),'Color','r','FontSize',8);
	end

	title({'One dimensional grid';' '},'FontSize',16,'FontWeight','bold');

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