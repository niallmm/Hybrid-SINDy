%% take averages across the instances  (3rd variable)
pathname = '';
mkdir(pathname);
changeplot
%%
hopyesavg = mean(hopyes, 3);
flyyesavg = mean(flyyes, 3);
condThetaavg = mean(condTheta,3);


%% Plots for paper


dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'SpringDynamics.fig'];
f = figure(14);
[C,h]=contourf(epsvec, clustersize,  hopyesavg');
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Model found for Spring Dynamics')
f.Children.XScale = 'log'; f.Children.YScale = 'log';
colormapnew
c = colorbar; 
% colorbar('off')
h.LevelStep = 0.2;
axis([1e-4 10 10 1.5e4])
f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
f.CurrentAxes.YTick = [1e1 1e2 1e3 1e4];

%%
savefig(plotname);
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'FlyDynamics.fig'];
f = figure(20);
[C,h]=contourf(epsvec, clustersize,  flyyesavg');
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Model found for Flying Dynamics')
f.Children.XScale = 'log'; f.Children.YScale = 'log';
colormapnew
c = colorbar; 
% colorbar('off')
h.LevelStep = 0.2;
axis([1e-4 10 10 1.5e4])
f.CurrentAxes.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
f.CurrentAxes.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondSpringDynamics.fig'];
f = figure(16);
contourf(epsvec, clustersize, log10(condThetaavg(:,:,:,1)'))
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Spring Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(condition #)';
hold on
caxis([0 6])
[C,h]=contour(epsvec, clustersize,  hopyesavg');
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondFlyDynamics.fig'];
f = figure(17);
contourf(epsvec, clustersize, log10(condThetaavg(:,:,:,2)'))
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Flying Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(condition #)';
hold on
caxis([0 6])
[C,h]=contour(epsvec, clustersize,  flyyesavg');
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
%%
savefig(plotname);
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('')
print([plotname(1:end-4) '.svg'],'-dsvg')
%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondEpsSpringDynamics.fig'];
f = figure(18);
contourf(epsvec, clustersize, log10(repmat(epsvec', [1 size(condThetaavg,2)])'.*condThetaavg(:,:,:,1)'))
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Spring Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(\epsilon \times condition #)';
hold on
caxis([-2 6])
[C,h]=contour(epsvec, clustersize,  hopyesavg');
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel(''); ylabel(''); title('');
print([plotname(1:end-4) '.svg'],'-dsvg')

%%
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'CondEpsFlyDynamics.fig'];
f = figure(19);
eps_cond =repmat(epsvec', [1 size(condThetaavg,2)])'.*condThetaavg(:,:,:,2)';
% contourf(epsvec, clustersize, log10(repmat(epsvec', [1 size(condThetaavg,2)])'.*condThetaavg(:,:,:,2)'))
contourf(epsvec, clustersize, eps_cond)
xlabel('noise, \epsilon')
ylabel('cluster size')
title('Flying Dynamics')
colormapnew
c = colorbar;
c.Label.String = 'k = log_{10}(\epsilon \times condition #)';
hold on
caxis([-2 6])
[C,h]=contour(epsvec, clustersize,  flyyesavg');
h.Parent.XScale = 'log'; h.Parent.YScale = 'log';
h.LineColor = [0.5 0.5 0.5];
h.LineStyle = '-';h.LineWidth = 1; h.LevelStep = 0.5;
axis([1e-4 10 10 1.5e4])
h.Parent.XTick = [1e-4 1e-3 1e-2 1e-1 1e0 1e1];
h.Parent.YTick = [1e1 1e2 1e3 1e4];
savefig(plotname);
%%
f.CurrentAxes.FontSize = 7.5;
f.CurrentAxes.LineWidth = 0.6;
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 4.8 2.5]*1.5;
f.PaperSize = [4.8 2.5]*1.5;
xlabel('')
ylabel('')
title('')
print([plotname(1:end-4) '.svg'],'-dsvg')

