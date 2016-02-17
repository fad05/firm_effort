function [vac, job, unemp, work] = valfun(wagefun, theta)
%Function calculates coefficients of linear value functions given some
%outside "wage function" (to calculates wage given any match quality pi) 
%and vector of market tightnesses "theta".
%vac, unemp: vectors of size Z (number of possible states), one number for
%every state.
%job,work: [Zx2] sized matrices, each row contains coefficients [a b] of
%function a*x+b.

% "truncated_normal_rule.m": Copyright 2014 John Burkardt, distributed under the GNU LGPL license.

%1)Pass global parameter structure.
global params

%1. Initial guess on coefficients.
v0 = 0.5*ones(params.nshocks,1);
u0 = v0;
j0 = [ones(params.nshocks,1) zeros(params.nshocks,1)];
w0 = j0;

%2.Value function iteration cycle.
delta = 10; %convergence criterion
while (delta>params.epsilon)
    %Given the value functions, find vector of threshold values of pi.
    pi_h_firm = (v0-j0(:,2))./j0(:,1);
    pi_h_worker = (u0-w0(:,2))./w0(:,1);
    %Calculate next iteration value functions
    %a. Calculate quadratures.
    
    vac = -params.kv + params.beta*((1-qmeet(theta)).*v0 + qmeet(theta));
    
end


%2. 



vac = feval(wagefun,3);
end