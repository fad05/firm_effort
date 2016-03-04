function f = qmeet(theta)
% "truncated_normal_rule.m": Copyright 2014 John Burkardt, distributed under the GNU LGPL license.
global params
f = max(0,min(1,params.acobbdoug*theta.^(-params.eta)));
end