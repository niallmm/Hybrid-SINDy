function dy = hopperflight(t,y)
% derivatives when hopper is in flight

dy = [ y(2,:);   % yp = velocity
       -1*ones(size(y(2,:)))]; % ypp = 1-kappa*y
    