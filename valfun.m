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
global params transmat

%1. Initial guess on coefficients.
v0 = 0.5*ones(params.nshocks,1);
u0 = v0;
j0 = [ones(params.nshocks,1) zeros(params.nshocks,1)];
w0 = j0;

%Initialize quadrature vectors;
quad_firm = zeros(params.nshocks,1);
quad_worker = zeros(params.nshocks,1);

% 2. Calculate quadrature points and weights.
[pts, wgt] = truncated_normal_rule(3,10,params.mu_pi,params.sigma_pi,params.a_pi, params.b_pi,'');

%3.Value function iteration cycle.
delta = 10; %convergence criterion
counter = 0;


while (delta>params.epsilon)
    counter = counter+1;
    %a. Calculate quadratures.
    for i = 1:params.nshocks
        func_firm = @(x) max(v0(i),j0(i,1)*x+j0(i,2));
        func_worker = @(x) max(u0(i), w0(i,1)*x+w0(i,2));
        quad_firm(i) = func_firm(pts)'*wgt;
        quad_worker(i) = func_worker(pts)'*wgt;
    end    
    %b. Calculate next iteration of value functions    
    %Vacancy VF
    vac = -params.kv + params.beta*transmat*((1-qmeet(theta)).*v0 + qmeet(theta).*quad_firm);
    %Unemployed VF
    unemp = params.beta*transmat*((1-pmeet(theta)).*u0 + pmeet(theta).*quad_worker);
    %Job VF
        % Part 1. Threshold values of pi.
        %In short: trying to find threshold values of pi. If they are
        %inside [0,1] interval, we check pihat and 1. If pihat is outside,
        %we check 0 and 1.        
        %Calculate threshold value of pi given the current guess on VF's.
        %We are going to check just two points: treshhold pi and pi = 1.
        %These two points fully determine a line.
        %To threshold pi's, identify zero and nonzero slopes of job and work valfuns
        pi_h_firm = zeros(params.nshocks,1); pi_h_worker = pi_h_firm;
        nonzero_j = (j0(:,1)~=0); zero_j = (j0(:,1)==0);
        nonzero_w = (w0(:,1)~=0); zero_w = (w0(:,1)==0);
        %For nonzero coefficients a.
        pi_h_firm(nonzero_j,1) = (v0(nonzero_j)-j0(nonzero_j,2))./j0(nonzero_j,1);
        pi_h_worker(nonzero_w,1) = (u0(nonzero_w)-w0(nonzero_w,2))./w0(nonzero_w,1);
        pi_h_firm(pi_h_firm>=1) = 0; %if pihat larger or equal than 1, just take it as 0
        pi_h_firm(pi_h_firm<0) = 0; %if pihat negative, just take it as 0
        % Part 2. Givet threshold values, compute next period J and W at
        % pihat(that can be 0) and 1.
        
        
        
        
    
end


%2. 



vac = feval(wagefun,3);
end