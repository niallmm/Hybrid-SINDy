function [value, isterminal, direction] = events_diverge(t,y);

value = max(abs(y))-200;
isterminal = 1;
direction = [];

