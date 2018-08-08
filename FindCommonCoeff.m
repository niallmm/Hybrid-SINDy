% Compare coefficients within a single model
% must have run FindMinModels first.
clear C1lib C2lib C1lib_ind C2lib_ind C1lib_hf C2lib_hf C1_hf_min C2_hf_min
clear C2_hf_x C1_hf_x C2_hf_y C1_hf_y
for freqnum = 1 % take the highest frequency model
    Modelstruct = HighFreq_II(:,freqnum); % get the structure of the highest frequency model
    coeff_tol = 1e-3; % tolerance of how far apart the coefficients can be
    
    % make a coefficient Library
    C1lib= uniquetol(C1', coeff_tol, 'ByRows', true)';
    C2lib= uniquetol(C2', coeff_tol, 'ByRows', true)';
    C1lib = C1lib(:,all(~isnan(C1lib),1)); % remove NaN
    C2lib = C2lib(:,all(~isnan(C2lib),1)); % remove NaN
    % find the ones that have the same structure as the highest frequency
    C1lib_ind = find(ismember(abs(C1lib')>0, Modelstruct(1:10)', 'rows'));
    C2lib_ind = find(ismember(abs(C2lib')>0, Modelstruct(11:20)', 'rows'));
    
    for kk = C1lib_ind'
        
        C1lib_hf = C1lib(:,kk)
        
        % find where these coefficients are
        C1hf_ind  = find(ismembertol(C1', C1lib_hf', coeff_tol, 'ByRows', true));
        C1_hf_min{kk,freqnum} = C1hf_ind'; % just store if these coefficients are active
        for ii = C1hf_ind'
            C1_hf_x(ii,kk) = Cluster{ii}.centroidx;
            C1_hf_y(ii,kk) = Cluster{ii}.centroidy;
        end
    end

    %%
    for kk = C2lib_ind'
        
        C2lib_hf = C2lib(:,kk)
        C2hf_ind  = find(ismembertol(C2', C2lib_hf', coeff_tol, 'ByRows', true));
         C2_hf_min{kk,freqnum} = C2hf_ind'; % just store if these coefficients are active
        for ii = C2hf_ind'
            C2_hf_x(ii,kk) = Cluster{ii}.centroidx;
            C2_hf_y(ii,kk) = Cluster{ii}.centroidy;
        end
    end

end
