function [] = plot(G,varargin)
% -----------------------------------------------------------------------------
% plot(varargin) plots the grid for debugging purposes.
% -----------------------------------------------------------------------------
	figure
	[hmin,hmax] = G.getGridSize();
	if size(varargin,2) == 0
		option.mode = 'debugging';
		option.NodeLabels = 'on';
		option.EdgeLabels = 'on';
		option.CellLabels = 'on';
	else
		option = varargin{1};
	end
	if strcmp(upper(option.mode),'DEBUGGING')
		hold on;
		numEdges = G.edges_.numEdges_;
		for k = 1:numEdges
			i1 = G.edges_.vList_(1,k);
			i2 = G.edges_.vList_(2,k);
			edges(:,2*(k-1)+1:2*k) = G.nodes_.xList_(:,G.edges_.vList_(1:2,k)');
			x = [G.nodes_.xList_(1,i1) G.nodes_.xList_(1,i2)];
			y = [G.nodes_.xList_(2,i1) G.nodes_.xList_(2,i2)];
			if G.edges_.sList_(2,k) > 0
				line(x,y,'Color','k');
			else
				line(x,y,'Color','r');
			end
			if strcmp(upper(option.EdgeLabels), 'ON')
				txt = sprintf('%d (%d)',k, G.edges_.hdgID_(k));
					text((x(1)+x(2))/2.0,(y(1)+y(2))/2.0,txt, ...
								'Color','k','FontSize',10,'FontWeight','bold');
			end
		end
		if strcmp(upper(option.NodeLabels),'ON')
			numNodes = G.nodes_.numNodes_;
			delta = 0.1;
			for k = 1:numNodes
				x = G.nodes_.xList_(1,k) + delta*hmin;
				y = G.nodes_.xList_(2,k) + delta*hmin;
				text(x,y,num2str(k),'Color','b','FontSize',10);
				plot(G.nodes_.xList_(1,k),G.nodes_.xList_(2,k),...
               'o','MarkerSize',10,...
						'MarkerEdgeColor','k','MarkerFaceColor','b');
			end
		end
		if strcmp(upper(option.CellLabels),'ON')
			numCells = G.cells_.numCells_;
			for k = 1:numCells
				gp = ( G.nodes_.xList_(:,G.cells_.vList_(1,k)) + ...
				       G.nodes_.xList_(:,G.cells_.vList_(2,k)) + ...
							 G.nodes_.xList_(:,G.cells_.vList_(3,k)))/3.0;
				text(gp(1),gp(2),num2str(k),'Color','r','FontSize',10);
			end
		end
	else
		hold on;
		cmap = colormap('jet');
		maxMatIndex = max(G.cells_(4,:));
		maxColorIndex = length(cmap);
		numCells = size(G.cells_,2);
		for k = 1:numCells
			Tk = G.cells_.vList_(1:3,k);
			colorIndex = G.cells_.mList_(k)*(maxColorIndex/maxMatIndex);
			patch(G.nodes_.xList_(1,Tk), G.nodes_.xList_(2,Tk),cmap(colorIndex,:),...
			      'EdgeColor','none');
		end
	end
	title({upper(G.file_);' '},'FontSize',11,'FontWeight','bold');
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