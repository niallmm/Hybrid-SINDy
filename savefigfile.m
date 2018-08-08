% code for saving figures of the right dimensions and font size for paper

% load plots
% savepath1 = 'C:\Users\niallmm\Dropbox\GitHub\SINDY_IC_paper\Figures\plotsforfigures\';
% mkdir(savepath1)
% open([path plotname]);

% change axis labels and range varies for each plot
% get from command line
if ~exist('xlimits')
xlimits = input('enter x limits [low high]');
ylimits = input('enter y limits [low high]');
xticks = input('enter ticks for x-axis [1 2 3]');
yticks = input('enter ticks for y-axis [1 2 3]');
end
axis([xlimits ylimits])
l = get(gca, 'Children');
xlabel('')
ylabel('')
title('')
ax  =gca;
% set(ax, 'FontSize', 5)
set(ax, 'FontSize', 7.5)
set(ax, 'FontName', 'Times New Roman')
if isprop(l, 'Marker')
    set(l,'MarkerSize',8)
    % set(l, 'MarkerSize', 2)
    set(l, 'Marker', '.')
end
set(l, 'linewidth', 0.6)
ax.XTick = xticks;
ax.YTick = yticks;
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';
ax.LineWidth = 0.8;
if exist('ylabelformat')
ax.YAxis.TickLabelFormat = ylabelformat;
end
fig = gcf;
fig.PaperUnits = 'centimeters';
ax.Units = 'centimeters';
% for 1/3 width figures
if exist('PP')
    fig.PaperPosition = PP;
    fig.PaperSize = PS;
else
    fig.PaperPosition = [0 0 3.5 2.5];
    fig.PaperSize = [3.5 2.5];
end
% % for 1/2 width figures 
% fig.PaperPosition = [0 0 8.5 6.5];
% fig.Papersize = [8.7 6.7];
% for full width figures
% fig.PaperPosition = [0 0 6 3.5];
% fig.PaperSize = [6.2 3.7];
% fig.PaperPosition = [0 0 3.25 2.25];
% fig.PaperSize = [3.25 2.25];

savefig(plotname);
print([plotname(1:end-4) '.svg'],'-dsvg')
% print(savename(1:end-4), '-dpdf')