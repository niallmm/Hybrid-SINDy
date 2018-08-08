function dy = hopperspring(t,y,kappa)
% derivatives when hopper is touching ground
 % yp = velocity
 % ypp = 1-kappa*y
dy = [ y(2,:);
    1- kappa.*(y(1,:)-1)]; 
    