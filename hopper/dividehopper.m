function [datafly, datahop] = dividehopper(d,dtol)
%% divide data into fly and springing
thresh = 1;
t_hop = d.tout(d.yout(:,1)<thresh-dtol);
y_hop = d.yout(d.yout(:,1)<thresh-dtol,:);
a_hop = d.aout(d.yout(:,1)<thresh-dtol,:);

t_fly = d.tout(d.yout(:,1)>thresh+dtol);
y_fly = d.yout(d.yout(:,1)>thresh+dtol,:);
a_fly = d.aout(d.yout(:,1)>thresh+dtol,:);

datafly.t = t_fly;
datafly.y = y_fly;
datafly.a = a_fly;
datahop.t = t_hop;
datahop.y = y_hop;
datahop.a = a_hop;