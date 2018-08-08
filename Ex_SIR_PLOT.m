pathname = 'C:\Users\niall\Dropbox\UWappointment\SINDy_Switching\SIR_plots\';
mkdir(pathname);

dateformatout = 'mmddyyyy';

%% Plot data for SIR model
% plots for paper figure
plotname = [pathname datestr(now, dateformatout) 'Svtime.fig'];
figure(116)
plot(tsave(1:365)/365, S_save(1:365), '.')
ylabel('susceptible')
xlimits = [0 1];
ylimits = [0 60];
xticks = [0 0.25 0.5 0.75 1];
yticks = [0 10 20 30 40 50 60];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile

figure(117)
plot(tsave(1:365)/365, I_save(1:365), '.')
ylabel('infected')
xlabel('time')
plotname = [pathname datestr(now, dateformatout) 'Ivstime.fig'];
savefigfile

%% Cluster Plot
close all
figure(15)
plot(x, y, 'o')
hold on
%
for i = 5:20:55 %1:numOfPts
        figure(15)
        plot(x(idx_xy(i,:)), y(idx_xy(i,:)),  '*')
        xlabel('S')
        ylabel('I')
        title('clustering in phase space')
        drawnow
end
plotname = [pathname datestr(now, dateformatout) 'datadrivencoordsSI.fig'];
xlimits = [0 60];
ylimits = [0 50];
xticks = 0:10:50;
yticks = 0:15:60;
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
savefigfile
%% 
 FindMinModels;
 
% plot the model frequency histogram
plotname = [pathname datestr(now, dateformatout) 'modelfreq.fig'];
xlimits = [0 40];%max(edges)];
ylimits = [0 max(counts)+20];
xticks = 0:5:80;
yticks = 0:200:max(counts)+20;
PP = [0 0 3.5 2.75]*0.8;
PS = [3.5 2.75]*0.8;
% PP = [0 0 4.5 3.5];
% PS = [0 0 4.5 3.5];
savefigfile

%% S vs I plot with found Coefficients overlayed
FindCommonCoeff;
close all
figure(15)
plot(x, y, 'o')
hold on

C1found_x = C1_hf_x(:,C1lib_ind(1))~=0;
plot(x(C1found_x), y(C1found_x),'x')
% plot(C1_hf_x(:,C1lib_ind(1)), C1_hf_y(:,C1lib_ind(1)),'*')
hold on
for ii = 2:length(C1lib_ind)
%     plot(C1_hf_x(:,C1lib_ind(ii)), C1_hf_y(:,C1lib_ind(ii)),'*')
    C1found_x = C1_hf_x(:,C1lib_ind(ii))~=0;
    plot(x(C1found_x), y(C1found_x),'x')
end
title('C1')
plotname = [pathname datestr(now, dateformatout) 'datadrivencoordsSI_C1.fig'];
xlimits = [0 60];
ylimits = [0 50];
xticks = 0:10:50;
yticks = 0:15:60;
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
%  savefigfile

figure(16)
plot(x, y, 'o')
hold on
C2found_x = C2_hf_x(:,C2lib_ind(1))~=0;
plot(x(C2found_x), y(C2found_x),'x')
% plot(C2_hf_x(:,C2lib_ind(1)), C2_hf_y(:,C2lib_ind(1)),'*')
hold on
for ii = 2:length(C2lib_ind)
    C2found_x = C2_hf_x(:,C2lib_ind(ii))~=0; 
    plot(x(C1found_x), y(C1found_x),'x')
end
title('C2')
plotname = [pathname datestr(now, dateformatout) 'datadrivencoordsSI_C2.fig'];
xlimits = [0 60];
ylimits = [0 50];
xticks = 0:10:50;
yticks = 0:15:60;
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
% savefigfile
%%
% plot coefficients in time
plotname = [pathname datestr(now, dateformatout) 'C1vtime.fig'];
figure(10)
semilogy(t_out/365,abs(C1)','*')
xlimits = [0 5];
ylimits = [1e-3 5];
xticks = 0:5;
yticks = [1e-3 1e-2 1e-1 1e0 1e1 1e2];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile

plotname = [pathname datestr(now, dateformatout) 'C2vtime.fig'];
figure(11)
semilogy(t_out/365,abs(C2)','*')
xxlimits = [0 5];
ylimits = [1e-3 1];
xticks = 0:5;
yticks = [1e-3 1e-2 1e-1 1e0 1e1 1e2];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile



%% plot aicmin vs time
plotname = [pathname datestr(now, dateformatout) 'minaicval.fig'];
figure
semilogy(HMt(1,:)/365, HMerror_minaic(1,:),'*')
 hold on
 for kk = 2:size(HMt, 1)
     plot(HMt(kk,:)/365, HMerror_minaic(kk,:),'*')
 end
xxlimits = [0 5];
ylimits = [5e-7 1];
xticks = 0:5;
yticks = [1e-6 1e-4 1e-2 1e0 1e2];
PP = [0 0 3.5 2.5]*0.8;
PS = [3.5 2.5]*0.8;
% savefigfile


%% plot time series
plotname = [pathname datestr(now, dateformatout) 'Svtime_model.fig'];

figure(116)
plot(tsave/365, S_save, 'o')
hold on
for freqnum = 1;
    for kk = C2lib_ind'
        C2ind = C2_hf_min{kk,freqnum};
        plot(tsave(C2ind)/365, S_save(C2ind),'*')
    end
end
%
ylabel('susceptible')
xlimits = [0 5];
ylimits = [0 60];
xticks = [0 1 2 3 4 5];
yticks = [0 10 20 30 40 50 60];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile
%
figure(117)
plot(tsave/365, I_save, 'o')
hold on
for freqnum = 1;
    for kk = C1lib_ind'
        C1ind = C1_hf_min{kk,freqnum};
        plot(tsave(C1ind)/365, I_save(C1ind),'*')
    end
end
%
plotname = [pathname datestr(now, dateformatout) 'Ivtime_model.fig'];
ylabel('susceptible')
xlimits = [0 5];
ylimits = [0 50];
xticks = [0 1 2 3 4 5];
yticks = [0 10 20 30 40 50];
PP = [0 0 3.5 1.5]*0.8;
PS = [3.5 1.5]*0.8;
savefigfile


%% relative AIC_c plot
% close all
modelnum = 1; %  this is the modelnumber in the HighFreq_II library
 % use this code for models 1 and 2, need to run for 3 and 4 each
 % before running next cell
plotname = [pathname datestr(now, dateformatout) 'AIC_rel1_' num2str(modelnum) '.fig'];
model1ind = find(min_aic_modelnum == modelnum);
for ii = model1ind
    AIC_rel = Cluster{ii}.aic_c - min(Cluster{ii}.aic_c);
    numcoeff = Cluster{ii}.numterms(1:length(AIC_rel));
     if any(numcoeff(find(0== AIC_rel)) == nnz(HighFreq_II(:,modelnum)))
        figure(2)
        semilogy(numcoeff, AIC_rel, 'ok')
     end
    hold on
end
%%
xlabel('number of terms')
ylabel('relative AICc')
xmax = 6;
axis([min(numcoeff) xmax 4 max(AIC_rel)])
patch([0 xmax xmax 0], ...
    [10 10  max(AIC_rel) max(AIC_rel)], ...
    [0.2 0.2 0.2], 'FaceAlpha', 0.6)
patch([0 xmax xmax 0], [4 4 7 7], [0.2 0.2 0.2], 'FaceAlpha', 0.3)
patch([0 xmax xmax 0], [0 0 2 2], [0.2 0.2 0.2], 'FaceAlpha', 0.1)
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
xticks = 0:1:xmax;
yticks = 1:4;
PP = [0 0 3.5 1]*0.8;
PS = [3.5 1]*0.8;
savefigfile
