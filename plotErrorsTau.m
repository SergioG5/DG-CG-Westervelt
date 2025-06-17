function plotErrorsTau(xv, yv, p)
% Plot errors in loglog scale (continuous line)
	nPlots =  p;
	marker = {'-s','o-','d-','p-','^-','*','s-','h-'};
	colors = {'r','b',[0,0.5,0],[0.5430 0 0],'magenta','cyan','r','magenta'};
	loglog(xv, yv, marker{nPlots}, 'Color', colors{nPlots}, 'MarkerFaceColor', colors{nPlots}, 'LineWidth', 2);
	for i=1:length(xv)-1
    	x = sqrt(xv(i)*xv(i+1));
		y = sqrt(yv(i)*yv(i+1));
		m = (log(yv(i+1))-log(yv(i)))/(log(xv(i+1))-log(xv(i)));
		text(x, y, sprintf('%.2f', m),'HorizontalAlignment','center','BackgroundColor','y','FontSize', 8, 'FontWeight', 'Bold');
	end
	box on
	grid on
	set(gca, 'LineWidth', 1.0, 'Fontsize', 12.0, 'FontWeight', 'bold');
	legend('q = 2', 'q = 3', 'q = 4', 'q = 5','Location','southeast');
	xlabel('$\log_{10}(N_{DoFs})$','Interpreter','LaTex', 'FontSize', 18);
	ylabel('');
end

