
function [xA, x0clust, tvectest] = buildinit_fromcluster(val_ntimes, ii, idx_xy,data, dt, k);
    [dlength, dnum] = size(data);
% dlength = length(x);
    % build a matrix that has the initial condition indices in first
    % comlumn and 1 in second, 2 in third and so on. up to the number
    % of timepoints we want
    indxM = idx_xy(ii,:)';
    tvectest= 0;
    for mm = 1:val_ntimes
        indxM = [indxM ones(k-1, 1)];
        tvectest = [tvectest mm*dt];
    end
    % use matrix multiplication to create a matrix where each row is
    % the indices for each initial condition and next val_ntimes points
    indxM2 = indxM*triu(ones(val_ntimes+1));
    % make sure you haven't gone past the end of the time series
    rowsM2 =find(all(indxM2<dlength,2));
    icused = length(rowsM2);
    % construc data matrix first index is x or y , 2nd is number of ...
    % time pts at each IC, 3rd is each IC
    clear MxA
    for kk = 1:dnum
%     MxA(1:dnum,:,:) = data(indxM2(rowsM2,:),:)';
        v = data(:,kk);
        MxA(kk,:,:) = v(indxM2(rowsM2,:))';
    end
%      MxA(kk,:,:) = x(indxM2(rowsM2,:))';
%      MxA(kk,:,:) = y(indxM2(rowsM2,:))';
    % put the data in the cluster and next time steps into a format
    % that the validateXi can read (treat each point in the cluster as
    % a different intial condition and stack x and xdt into a column
    % vector). So each initial condition lives in a cell.
    xA = num2cell(MxA, [2 1]);
    x0clust = data(indxM2(rowsM2,1),:)';
%         x0clust =  [ x(indxM2(rowsM2,1))'; y(indxM2(rowsM2,1))'];