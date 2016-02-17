%Create parameter structure

% "truncated_normal_rule.m": Copyright 2014 John Burkardt, distributed under the GNU LGPL license.

global params
params = struct('nshocks', 3, ...
'beta', 0.9966,...
'yb', 1,...
'yg', 1.9, ...
'adisutil',1.0, ...
'lambda', 0.0085, ... %exogenous destruction
'mu_pi', 0.0, ... %mean of parent normal distribution
'sigma_pi', 0.32, ... %SD of parent normal distribution
'a_pi', 0.0, ... %left truncation of parent normal distribution
'b_pi', 1.0, ... %right truncation of parent normal distribution
'zeta', 0.5, ...
'eta', 0.5, ...
'kx', 2, ...
'kv', 0.5566, ...
'acobbdoug', 0.6, ...
'epsilon', 1E-4, ...
'alpha', 0.13,...
'nquad', 10); %number of quadrature points
