function plotErrorsP(xv, yv, p)
% Plot errors in loglog scale (continuous line)
	nPlots =  p;
	marker = {'-s','o-','d-','p-','^-','*','s-','h-'};
	colors = {'r','b',[0,0.5,0],[0.5430 0 0],'magenta','cyan','r','magenta'};
	semilogy(xv, yv, marker{nPlots}, 'Color', colors{nPlots}, 'MarkerFaceColor', colors{nPlots}, 'LineWidth', 2);
	box on
	grid on
	set(gca, 'LineWidth', 1.0, 'Fontsize', 12.0, 'FontWeight', 'bold');
end

