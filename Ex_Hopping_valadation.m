% validation for hopping model
% run after Ex_Hopping.m

% =========================================================================
% run new data for validation step.
% =========================================================================
options_test = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events', @events_diverge);

% start with mass slightly below it's equilibrium position
% velocity downward
yinitvec_val = [0.84 -0.11
    0.77 -0.12
    0.83 -0.13
    0.79 -0.10
    0.82 -0.11];
p.yinitvec = yinitvec_val;

[data]  = RUN_Hopper(p);
t_val = data.tout; 
y_val = data.yout+ rand(size(data.yout))*p.eps; 
a_val =data.aout;


%% Find cluster in validation data that is nearest the centroid of the original clusters
% initialize a structure library
II_Xi_aic =[]; ind_all_aic = [];

for ii = 1:length(Cluster)
    Xicomb = Cluster{ii}.Xistruct.Xicomb;
    clear x2 y2
    x2= y_val(1:end-1,1); % height
    y2= y_val(1:end-1,2); % velocity
    
    dimMan = 2; % dimension of the Manifold
    numOfPts = length(x2); % number of points in the data
    k = 30; % number of points in a cluster
    
    idx_xy2 = knnsearch([x2 y2],...
        [Cluster{ii}.centroidx Cluster{ii}.centroidy],'K',k);
    idx_xy2 = idx_xy2(:,2:end);

    val_ntimes = floor(1.5/dt); %length of comparision not including ICs
    data = [x2 y2];
    [xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes, 1, idx_xy2, data, dt, k);
    
    val.tA =tvectest;
    val.xA =xA;
    val.x0 =x0clust;
    val.options = options_test;
    
    clear IC aic_c aic_cv tdivmax tdivmin nodynam_any tdivavg tdiv 
    clear nodynam abserror_bf_switch mean_e numterms
    for kk = 1:length(Xicomb)
        disp(['Cluster num:' num2str(ii) ' of' num2str(length(Cluster))])
        disp(['Model num:' num2str(kk) ' of ' num2str(length(Xicomb)) ' for this cluster'])
        Xi = Xicomb{kk}
        
        clear  savetB savexB
        [savetB, savexB] = validation_sims(Xi, Thetalib,val, plottag);
        
        for ll = 1:length(xA)
            clear tdiv1 nodynam1 abserror abserror_avg1 RMSE xAcomp xBcomp
            [tlength, nterms] = size(savexB{ll});
            [tdiv1, abserror, abserror_avg1, RMSE] = calc_tlength(xA{ll}, savexB{ll},val);
            tdiv(ll,:) = tdiv1;
            abserror_bf_switch(ll,:) = abserror_avg1;
        end
        abserror = [abserror_bf_switch(:,2); abserror_bf_switch(:,1)];
        IC{kk} = ICcalculations(abserror, nnz(Xi),size(x0clust,2));
         aic_c(kk,:) = IC{kk}.aic_c;
         numterms(kk) = nnz(Xi);
        tdivmax(kk,:) = max(tdiv);
        tdivmin(kk,:) = min(tdiv);
        tdivavg(kk,:) = mean(tdiv);
        mean_e(kk,:) = mean(abserror_bf_switch);
        valsims{kk}.savetB = savetB;
        valsims{kk}.savexB = savexB;
        valsims{kk}.val = val;
        
    end

    minIC = min(aic_c);
    Xi_select = find( abs(aic_c-minIC)<3)
    Xicomb{Xi_select}
    % Make a library of models below AIC_c threshold
    [II_Xitemp, II_ind] = build_Xi_Library(Xicomb(Xi_select), II_Xi_aic);
    II_Xi_aic = II_Xitemp;
    ind_all_aic = [ind_all_aic; II_ind]
    
    Cluster{ii}.tdivmax = tdivmax;
    Cluster{ii}.abserror_bf = mean_e;
    Cluster{ii}.Xi_ind = Xi_select;
    Cluster{ii}.IC = IC;
    Cluster{ii}.aic_c =aic_c;
    Cluster{ii}.numterms =numterms;
    Cluster{ii}.valsims = valsims;
    
end

dateformatout = 'mmddyyyy';
save([datestr(now, dateformatout) 'Hopping_validation_eps1e_6.mat'])

