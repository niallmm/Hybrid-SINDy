
% Find high frequency model numbers
edges = 0.5:1:length(II_Xi_aic)+0.5;
figure(25)
h = histogram(ind_all_aic, edges);
title('Number of clusters with each model')
xlabel('model #')

counts = h.Values;
[B, I] = sort(counts,'descend'); % sort the models by their count

HighFreq_II = []; % make library of high freqency models
for kk = 1:6
    I(kk) % index of the kk frequent model
    modelfreq = B(kk)
    HighFreq_II = [HighFreq_II II_Xi_aic(:,I(kk))];
end

%% for each cluster
jj = [];
for ii = 1:length(Cluster)
    
    % check that there are models and aic_c scores for this cluster
    if ~isfield(Cluster{ii},'Xi_ind')||~isfield(Cluster{ii},'aic_c')
        disp(['no data for Cluster(' num2str(ii) ').Xi_ind'])
        model_ind(ii) = NaN;
    else
        % see if any of the models in this cluster are in our High Freq
        % library
        Xicomb = Cluster{ii}.Xistruct.Xicomb; % model library for this cluster
        Xiaic_loc = Xicomb(Cluster{ii}.Xi_ind); % library of just the models with low aic
        clear Xiind inLib XiMat
        for kk = 1:length(Xiaic_loc) % for each unique returned coefficient matrix
            XiMataic(:,kk) = reshape(Xiaic_loc{kk}, [],1);
        end
        IItestaic = unique((abs(XiMataic)>0)','rows', 'stable')'; % make a structure vector of 0s and 1s
        inLib = find(ismember(IItestaic', HighFreq_II', 'rows')); 
        %indices of low aic models that are in HighFreq_II library
        
        %want to find the index of the low aic, high frequency models for
        %the models this cluster (Xicomb) 
        for kk = 1:length(Xicomb)
            XiMat(:,kk) = reshape(Xicomb{kk},[],1);
        end
        IItest = unique((abs(XiMat)>0)','rows','stable')'; % make a structure vector of 0s and 1s

        
        for kk=inLib'
            Xiind = find(ismember(IItest', IItestaic(:,kk)', 'rows'));
            Xiind % index in the model library for this cluster
            jj = find(ismember(HighFreq_II', IItestaic(:,kk)','rows'));
            % jj is the index of the current model in the highfrequency
            % library
            
            if length(Xiind)>0 % make sure there is a model to store info on
                Xicomb{Xiind}
                Cluster{ii}.aic_c(Xiind,:)
                HMaic_c(jj,ii) = mean(Cluster{ii}.aic_c(Xiind, :));
                HMerror(jj,ii) = mean(Cluster{ii}.abserror_bf(Xiind, :));
                HMtdivmax(jj,ii) = min(Cluster{ii}.tdivmax(Xiind, :));
                HMx(jj,ii) = Cluster{ii}.centroidx;
                HMy(jj,ii) = Cluster{ii}.centroidy;
                if exist('t_out')
                    HMt(jj,ii) = t_out(ii);
                end
            end
            
        end
        % if there are low rel(AIC_c) models that are also in the
        % high-frequency library then find the minimum AIC_c one
        if ~isempty(jj)&&length(HMaic_c) >=ii % check to make sure they exist
            if any(HMaic_c(:,ii)) % if at least one of the high frequency models has an aic_c calculated
                min_aic_modelnum(ii) = find(min(HMaic_c(abs(HMaic_c(:,ii))>0,ii)) == HMaic_c(:,ii));
                HMxmin(min_aic_modelnum(ii),ii) = Cluster{ii}.centroidx;
                HMymin(min_aic_modelnum(ii),ii) = Cluster{ii}.centroidy;
                HMerror_minaic(min_aic_modelnum(ii),ii) = HMerror(min_aic_modelnum(ii),ii);
                min_aic_Xi_ind = find(ismember(IItest', HighFreq_II(:,min_aic_modelnum(ii))', 'rows'));
                C1(:,ii) = Xicomb{min_aic_Xi_ind}(:,1);
                C2(:,ii) = Xicomb{min_aic_Xi_ind}(:,2);
               if exist('t_out')
                    HMt_min(min_aic_modelnum(ii),ii) = t_out(ii);
               end
            else
                min_aic_modelnum(ii) =NaN;
            end
        else
            min_aic_modelnum(ii) = NaN;
            HMxmin(:,ii) = NaN;
            HMymin(:,ii) = NaN;
            C1(:,ii) =NaN*ones(length(Xicomb{1}),1);
            C2(:,ii) =NaN*ones(length(Xicomb{1}),1);
        end
    end
    
end
