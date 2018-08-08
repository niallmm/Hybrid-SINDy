function [tdiv, abserror, abserror_avg, RMSE] = calc_tlength(xA, xB, val)

        % define vectors to compare, making sure they are the same
        % length, and calculate length
        [tlength, nterms] = size(xB);
        xAcomp =xA(1:tlength, :);
        xBcomp = xB;
        
        % calculate errors
        if length(xAcomp) == length(xB)
            error1 = xAcomp-xBcomp;
        else
            error('model and validation time-series are not the same length')
        end
        
        abserror = abs(error1);
        tvec = val.tA;
        
        for ll=1:nterms % for each term
            % find the time at which the error changes significantly
            % (switching point)
            if size(error1,1)>1
                switch_ind = findchangepts(error1(:,ll), 'Statistic','rms');
            else
                switch_ind = [];
            end
            
            if any(size(switch_ind)<1)% check if it's empty
                if any(size(switch_ind) == 0)
                 switch_ind = length(abserror(:,ll));
                end
            end
            tdiv(ll) = tvec(switch_ind);

            if length(xAcomp) == length(xBcomp)
                abserror_avg(ll) = sum(abs(xAcomp(1:switch_ind,ll)-xBcomp(1:switch_ind,ll)))/switch_ind;
                RMSE(ll) = sqrt(sum(abs(xAcomp(1:switch_ind,ll)-xBcomp(1:switch_ind,ll)).^2)/switch_ind);
            else
                error('model and validation time-series are not the same length')
            end
        end
        
        
       