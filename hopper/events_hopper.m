function [value, isterminal, direction] = events_hopper(t,y)
% y(2) is velocity
% y(1) is position
value = y(1)-1;
isterminal = 1;
direction = 0;