function f = pmeet(theta)
% "truncated_normal_rule.m": Copyright 2014 John Burkardt, distributed under the GNU LGPL license.
global params
f = max(0,min(1,params.acobbdoug*theta.^(1-params.eta)));
end