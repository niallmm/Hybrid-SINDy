function [savetB, savexB] = validation_sims(Xi, Thetalib,validation, plottag)
% perform a cross validation given Xi
% returns time vector savetB and validation series data savexB
% these structures contain time-series data for each

polyorder = Thetalib.polyorder;
usesine = Thetalib.usesine;

x0 = validation.x0;
tA = validation.tA;
xA = validation.xA;
options = validation.options;
% null1= validation.null1;
% null2 = validation.null2;

if plottag ==2
    Xi
end

[nterms, ninits] = size(x0);
for kk = 1:ninits % for each initial condition
    if iscell(tA)
        tvec = tA{kk};
    elseif isnumeric(tA)
        tvec = tA;
    end
    
    sol=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),...
        tvec,x0(:,kk),options);  % approximate
    % if the solution diverges quickly:
    dydx = abs(sol.y(end)-sol.y(end-1)/(sol.x(end)-sol.x(end-1)));
    
    
    % check if the test model produced a time series longer than
    % the validation time-series or vice-versa
    ind_tlast = find(tvec<=(max(sol.x)), 1, 'last');
    ind_tlast2 = find(sol.x<=(max(tvec)), 1, 'last');
    
    if (abs(sol.y(end))>=1e5)||(dydx>1e6); %if the last point is larger than divergence criteria
        savetB{kk} = tvec;
        savexB{kk} = [sol.y(:,1)'; sol.y(end)*ones(length(tvec)-1,nterms)];
        ind_tlast = length(tvec);
    elseif tvec(ind_tlast)<sol.x(ind_tlast2) % if the validation time-series
        % is longer than the test model, evaluate the test model at
        % the validation time-series points, until the last time
        % point that is calculated by both
        savetB{kk} = tvec(1:ind_tlast);
        savexB{kk} = deval(sol, tvec(1:ind_tlast))';
        
    elseif tvec(ind_tlast)>=sol.x(ind_tlast2) % if the test model results are
        % longer than or equal to the validation time-series,
        % evaluate at validation points only
        savetB{kk} = tvec;
        savexB{kk} = deval(sol, tvec)';
    else
        disp('some weird case!!!')
    end
    
    
    % plot the validation data and current model
    if (plottag ==1)
        Xi;
        figure(111)
        plot(tvec(1:ind_tlast), xA{kk}(1:ind_tlast,:), 'o')
        hold on
        plot(tvec(1:ind_tlast), savexB{kk}(1:ind_tlast,:))
        xlabel('time')
        ylabel('state variables')
        title('validation')
        drawnow
        hold off
    end
    
    
    
end



