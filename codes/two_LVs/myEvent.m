function [value, isterminal, direction] = myEvent(t, y)

value      = (y(3) < 0.999);
isterminal = 1;   % Stop the integration
direction  = 0;