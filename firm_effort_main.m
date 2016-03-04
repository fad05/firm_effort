% Created by Alex Filatov 16.02.2016
% Main file to the "Firm Effort" project
% "truncated_normal_rule.m": Copyright 2014 John Burkardt, distributed under the GNU LGPL license.
clear all
global params transmat zprod
parameters
%For now, random transition matrix
% transmat = rand(params.nshocks);
% transmat = transmat./repmat(sum(transmat,2),1,params.nshocks);
transmat = [1,0,0;0,1,0;0,0,1];
%Some matrix of states
%zprod = linspace(0.9,1.1,params.nshocks)';
zprod = [1;1;1];
theta = [1; 1; 1];
%a = 1+rand(params.nshocks,1);
a = 0.4432*ones(params.nshocks,1);
b = 1.3135*ones(params.nshocks,1);
wagefun = @(x) wage(x,a,b);
c = valfun(wagefun,theta);

