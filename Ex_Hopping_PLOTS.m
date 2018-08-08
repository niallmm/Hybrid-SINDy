% input a save location in the line below:
pathname = '';
mkdir(pathname);

dateformatout = 'mmddyyyy';

%% plots for paper figure Hopper
plotname = [pathname datestr(now, dateformatout) 'positionvtime.fig'];
figure(116)
plot(t_out,y_out(:,1),'.')
ylabel('position')
xlimits = [0 6];
ylimits = [0.75 1.5];
xticks = [0 1 2 3 4 5 6];
yticks = [0.5 0.75 1 1.25 1.5];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile

figure(117)
plot(t_out,y_out(:,2),'.')
ylabel('velocity')
xlabel('time')
plotname = [pathname datestr(now, dateformatout) 'velocityvtime.fig'];
xlimits = [0 6];
ylimits = [-1 1];
xticks = [0 1 2 3 4 5 6];
yticks = [-1 -0.5 0 0.5 1];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile

%% Cluster Plots
close all
figure(15)
plot(x, y, 'o')
hold on

for i = 5:10:25 %numOfPts
%     if plottag >2
        figure(15)
        plot(x(idx_xy(i,:)), y(idx_xy(i,:)), '*')        
        xlabel('position')
        ylabel('velocity')
        title('clustering in phase space')
        drawnow
%     end
end
dateformatout = 'mmddyyyy';
plotname = [pathname datestr(now, dateformatout) 'datadrivencoords.fig'];
xlimits = [0.6 1.6];
ylimits = [-1 1];
xticks = [0.6 0.8 1 1.2 1.4 1.6];
yticks = [-1 -0.5 0 0.5 1];
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
savefigfile

%%
% plot some of the coefficients
for clustnum = [5 14 25]
    delta_aic_c= Cluster{clustnum}.aic_c-min(Cluster{clustnum}.aic_c)
    nLib = length(Cluster{clustnum}.Xistruct.Xicomb);
    figure(clustnum)
    ii =1;
    dateformatout = 'mmddyyyy';
    plotname = [pathname datestr(now, dateformatout) 'Coeff' num2str(clustnum) '.fig']
    for  kk = 1:nLib
        
        vc = Cluster{clustnum}.Xistruct.Xicomb{kk}(:,1)
        ac = Cluster{clustnum}.Xistruct.Xicomb{kk}(:,2)
        
        
        subplot(nLib,2,ii)
        b1 = bar(abs(vc));
        box off
        b1.BarWidth = 0.1;
        b1.Parent.YScale = 'log';
        b1.Parent.YColor = [1 1 1];
        b1.Parent.XAxis.TickLength = [0 0];
        b1.Parent.XAxis.Limits = [0.5 6.5];
        b1.Parent.XAxis.LineWidth=1;
        b1.Parent.FontSize = 6;
        b1.Parent.Parent.PaperPosition = [0 0 3.5 2.5];
        b1.Parent.Parent.PaperSize = [3.5 2.5];
        if  kk ~= nLib
            b1.Parent.XAxis.TickValues = [];
        end
        ii = ii+1;
        subplot(nLib,2,ii)
        b2 = bar(abs(ac));
        b2.BarWidth = 0.1;
        box off
        b2.Parent.YScale = 'log';
        b2.Parent.YColor = [1 1 1];
        b2.Parent.XAxis.TickLength = [0 0];
        b2.Parent.XAxis.Limits = [0.5 6.5];
        b2.Parent.XAxis.LineWidth=1;
        b2.Parent.FontSize = 6;
        b2.Parent.Parent.PaperPosition = [0 0 3.5 2.5];
        b2.Parent.Parent.PaperSize = [3.5 2.5];
        if  kk ~= nLib
            b2.Parent.XAxis.TickValues = [];
        end
        ii= ii+1;
        
    end
    savefig(plotname);
    print([plotname(1:end-4) '.svg'],'-dsvg')
end


%% Figure of validation data clustering

figure(119)
plot(y_val(:,1), y_val(:,2),'o')
hold on
plot(y_out(:,1), y_out(:,2),'*')

plotname = [pathname datestr(now, dateformatout) 'validation_datadrivencoords.fig'];
close all
figure(15)
plot(y_val(:,1), y_val(:,2),'o')
hold on
xlabel('position')
ylabel('velocity')
title('validation')

for ii = [5 15 25]
    idx_xy2 = knnsearch([x2 y2],...
        [Cluster{ii}.centroidx Cluster{ii}.centroidy],'K',k);
    idx_xy2 = idx_xy2(:,2:end);
    figure(15)
    plot(x2(idx_xy2), y2(idx_xy2), '*')
    plot(Cluster{ii}.centroidx, Cluster{ii}.centroidy, 'x')
end

xlimits = [0.6 1.6];
ylimits = [-1 1];
xticks = [0.6 0.8 1 1.2 1.4 1.6];
yticks = [-1 -0.5 0 0.5 1];
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
savefigfile

%% Figure of validation time series (will plot whatever the latest validated cluster was)
jj = 5;
valsims = Cluster{jj}.valsims;
Xicomb = Cluster{jj}.Xistruct.Xicomb;
for kk = 9 % correct model for Cluster 4
    
    Xi = Xicomb{kk}
    clear  savetB savexB
    savetB = valsims{kk}.savetB;
    savexB = valsims{kk}.savexB;
    idx_xy2 = knnsearch([x2 y2],...
        [Cluster{jj}.centroidx Cluster{jj}.centroidy],'K',k);
    idx_xy2 = idx_xy2(:,2:end);
    val_ntimes = floor(1.5/dt); %length of comparision not including ICs
    data = [x2 y2];
    [xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes, 1, idx_xy2, data, dt, k);
%     [savetB, savexB] = validation_sims(Xi, Thetalib,val, plottag);
    for ii = 1:size(xA, 3)
        if ii == 1
            figure(115)
            plot(savetB{ii}, savexB{ii}(:,1),'-c')
            hold on
            plot(tvectest,xA{ii}(:,1),'.b')
            ylabel('position')
            drawnow
        end
        
        figure(116)
        plot(savetB{ii}, savexB{ii}(:,1),'-c')
        hold on
        plot(tvectest,xA{ii}(:,1),'.b')
        ylabel('position')
        drawnow
        
        figure(117)
        plot(savetB{ii}, savexB{ii}(:,2),'-c')
        hold on
        plot(tvectest,xA{ii}(:,2),'.b')
        ylabel('velocity')
        xlabel('time')
        drawnow
        
    end
    
    
end
%%
plotname = [pathname datestr(now, dateformatout) 'positionvtime_val1.fig'];
figure(115)
xlimits = [0 1.5];
ylimits = [0.75 1.6];
xticks = [0 0.5 1  1.5];
yticks = [0.75 1 1.25 1.5];
PP = [0 0 3.5 2]*0.8;
PS = [3.5 2]*0.8;
savefigfile
%%
plotname = [pathname datestr(now, dateformatout) 'positionvtime_val.fig'];
figure(116)
xlimits = [0 1.5];
ylimits = [0.75 1.6];
xticks = [0 0.5 1 1.5];
yticks = [0.5 1 1.5];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile

plotname = [pathname datestr(now, dateformatout) 'velocityvtime_val.fig'];
figure(117)
xlimits = [0 1.5];
ylimits = [-0.5 1];
xticks = [0 0.5 1 1.5];
yticks = [-0.5 0 0.5 1];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile
%% relative AIC_c plot
close all
modelnum = 4; %  this is the modelnumber in the HighFreq_II library
 % use this code for models 1 and 2, need to run for 3 and 4 each
 % before running next cell
plotname = [pathname datestr(now, dateformatout) 'AIC_rel1_' num2str(modelnum) '.fig'];
% model1ind = find(min_aic_modelnum == modelnum);
for ii = [5 14 25]
    AIC_rel = Cluster{ii}.aic_c - min(Cluster{ii}.aic_c);
    numcoeff = Cluster{ii}.numterms(1:length(AIC_rel));
        figure(ii)
        semilogy(numcoeff, AIC_rel, 'ok')
    hold on
end
%%
xlabel('number of terms')
ylabel('relative AICc')
xmax = 14;
axis([min(numcoeff) xmax 4 max(AIC_rel)])
% patch([0 xmax xmax 0], ...
%     [10 10  max(AIC_rel) max(AIC_rel)], ...
%     [0.2 0.2 0.2], 'FaceAlpha', 0.6)
% patch([0 xmax xmax 0], [4 4 7 7], [0.2 0.2 0.2], 'FaceAlpha', 0.3)
% patch([0 xmax xmax 0], [0 0 2 2], [0.2 0.2 0.2], 'FaceAlpha', 0.1)
%
xlimits = [0 xmax];
ylimits = [4 max(AIC_rel)];
xticks = 0:2:xmax;
yticks = [5:5:15 30 50:50:300];
PP = [0 0 3.5 1.75]*0.8;
PS = [3.5 1.75]*0.8;
savefigfile
%% run previous cell first
modelnum = 1; %  this is the modelnumber in the HighFreq_II library
 % use this code for models 1 and 2, need to run for 3 and 4 each
 % before running next cell
 model1ind = find(min_aic_modelnum == modelnum);
plotname = [pathname datestr(now, dateformatout) 'AIC_rel2_' num2str(modelnum) '.fig'];
figure(3);
for ii = model1ind
    AIC_rel = Cluster{ii}.aic_c - min(Cluster{ii}.aic_c);
    numcoeff = Cluster{ii}.numterms(1:length(AIC_rel));
    if any(numcoeff(find(0== AIC_rel)) == nnz(HighFreq_II(:,modelnum)))
        plot(numcoeff, AIC_rel, 'or')
        hold on
    end
end
%%
axis([min(numcoeff) xmax 0 4])
%    axis([min(numcoeff) 20 0 max(AIC_rel)])
ylabel('relative AICc')
patch([0 xmax xmax 0], ...
    [10 10  max(AIC_rel) max(AIC_rel)], ...
    [0.2 0.2 0.2], 'FaceAlpha', 0.6)
patch([0 xmax xmax 0], [4 4 7 7], [0.2 0.2 0.2], 'FaceAlpha', 0.3)
patch([0 xmax xmax 0], [0 0 2 2], [0.2 0.2 0.2], 'FaceAlpha', 0.1)

xlimits = [0 xmax];
ylimits = [0 4];
xticks = 0:2:xmax;
yticks = 1:4;
PP = [0 0 3.5 1]*0.8;
PS = [3.5 1]*0.8;
savefigfile

%% Run script to find Min Models

FindMinModels;
% plot the model frequency histogram
plotname = [pathname datestr(now, dateformatout) 'modelfreq.fig'];
xlimits = [0 max(edges)];
ylimits = [0 max(counts)+2];
xticks = 0:5:30;
yticks = 0:50:max(counts)+2;
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
% PP = [0 0 4.5 3.5];
% PS = [0 0 4.5 3.5];
savefigfile
%% plot of lowest aic model for each cluster
figure
plot(HMxmin(1,:), HMymin(1,:),'*')
hold on
plot(HMxmin(2,:), HMymin(2,:),'*')
plot(HMxmin(3,:), HMymin(3,:),'*')
hold on
plot(HMxmin(4,:), HMymin(4,:),'*')

%
axis([0.6 1.6 -1 1])
plotname = [pathname datestr(now, dateformatout) 'model_phasespace.fig'];
xlimits = [0.6 1.6];
ylimits = [-1 1];
xticks = [0.6 0.8 1 1.2 1.4 1.6];
yticks = [-1 -0.5 0 0.5 1];
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
% PP = [0 0 4.5 3.5];
% PS = [4.5 3.5];
savefigfile

%%
% plot coefficients
C2nz = C2; C2nz(C2==0) = NaN;
figure
plot(t_out(1:end-1),C2nz','*')
xlimits = [0 6];
ylimits = [min(min(C2))  max(max(C2)) ];
yticks = round(min(min(C2)),-1):100:round(max(max(C2)),-1)+50
xticks = [0 1 2 3 4 5 6];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
plotname = [pathname datestr(now, dateformatout) 'C2vtime.fig'];
savefigfile
    
     ylimits = [-11 11]
     yticks = -12:2:12
     plotname = [pathname datestr(now, dateformatout) 'C2_' num2str(ii) 'vtimezoom.fig'];
     savefigfile  
     %%
%      
% for ii =  1:size(C2,1)
%     plotname = [pathname datestr(now, dateformatout) 'C2_' num2str(ii) 'vtime.fig'];
% 
%     % plot(t_out(1:end-1),C2','*')
% %     if any(ii == [1 2])
% %         plot(t_out(1:end-1),abs(C2(ii,:))','*')
% %     else
%         plot(t_out(1:end-1),C2(ii,:)','*')
%         yticks = round(min(C2(ii,:)),-2):100:round(max(C2(ii,:)),-2)
% %     end
%     xlimits = [0 6];
%     ylimits = [min([-10 C2(ii,:)])  max([C2(ii,:) 10] )];
%     xticks = [0 1 2 3 4 5 6];
%     
% %     yticks = [5e0 1e1 1e2];
%     PP = [0 0 3.5 1]*0.8;
%     PS = [3.5 1]*0.8;
% %     savefigfile
% end
%% plot aicmin vs time
t1ind = length(HMt_min)/3;
plotname = [pathname datestr(now, dateformatout) 'minaicval.fig'];
figure
semilogy(HMt(1,1:t1ind), HMerror_minaic(1,1:t1ind),'*')
hold on
plot(HMt(2,1:t1ind), HMerror_minaic(2,1:t1ind),'*')
plot(HMt(3,1:t1ind), HMerror_minaic(3,1:t1ind),'*')
plot(HMt(4,1:t1ind), HMerror_minaic(4,1:t1ind),'*')
axis([0 6  0.7*min(min(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0))) 1.5*max(max(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0)))])
xlimits = [0 6];
ylimits = [0.7*min(min(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0))) 1.5*max(max(HMerror_minaic(HMerror_minaic(:,1:t1ind)>0)))];
xticks = [0 1 2 3 4 5 6];
yticks = [1e-6 1e-4 1e-2 1e0 1e2 1e4];
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
savefigfile


%% plot time series
t1ind = length(HMt_min)/3;

plotname = [pathname datestr(now, dateformatout) 'positionvtime_model.fig'];
figure(116)
plot(HMt_min(1,1:t1ind), HMxmin(1,1:t1ind), '.')
hold on
plot(HMt_min(2,1:t1ind), HMxmin(2,1:t1ind), '.')
plot(HMt_min(3,1:t1ind), HMxmin(3,1:t1ind), '.')
plot(HMt_min(4,1:t1ind), HMxmin(4,1:t1ind), '.')
ylabel('position')
xlimits = [0 6];
ylimits = [0.75 1.5];
xticks = [0 1 2 3 4 5 6];
yticks = [0.5 0.75 1 1.25 1.5];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile

figure(117)
plot(HMt_min(1,1:t1ind), HMymin(1,1:t1ind), '.')
hold on
plot(HMt_min(2,1:t1ind), HMymin(2,1:t1ind), '.')
plot(HMt_min(3,1:t1ind), HMymin(3,1:t1ind), '.')
plot(HMt_min(4,1:t1ind), HMymin(4,1:t1ind), '.')
ylabel('velocity')
xlabel('time')
plotname = [pathname datestr(now, dateformatout) 'velocityvtime_model.fig'];
xlimits = [0 6];
ylimits = [-1 1];
xticks = [0 1 2 3 4 5 6];
yticks = [-1 -0.5 0 0.5 1];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile


