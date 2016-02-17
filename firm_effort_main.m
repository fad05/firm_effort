% Created by Alex Filatov 16.02.2016
% Main file to the "Firm Effort" project
% "truncated_normal_rule.m": Copyright 2014 John Burkardt, distributed under the GNU LGPL license.
clear all
global params transmat zprod
parameters
%For now, random transition matrix
transmat = rand(params.nshocks);
transmat = transmat./repmat(sum(transmat,2),1,params.nshocks);
%Some matrix of states
zprod = linspace(0.9,1.1,params.nshocks)';
theta = [0.95; 1; 1.05];
a = 1+rand(params.nshocks,1);
b = zeros(params.nshocks,1);
wagefun = @(x) wage(x,a,b);
c = valfun(wagefun,theta);

