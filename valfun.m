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
global params transmat zprod

%1. Initial guess on coefficients.
v0 = 0.5*ones(params.nshocks,1);
u0 = 0.1*v0;
j0 = [ones(params.nshocks,1) zeros(params.nshocks,1)];
w0 = j0;

%Initialize quadrature vectors;
quad_firm = zeros(params.nshocks,1);
quad_worker = zeros(params.nshocks,1);

% 2. Calculate quadrature points and weights.
[pts,wgt] = tauchen(params.a_pi,params.b_pi,params.sigma_pi,params.nquad);

%3.Value function iteration cycle.
delta = 10; %convergence criterion
counter = 0;

tic
while (delta>params.epsilon)
    counter = counter+1;
    %a. Calculate quadratures.
    jhandle = @(x)max(repmat(v0,1,params.nquad),j0(:,1)*x'+repmat(j0(:,2),1,params.nquad));
    whandle = @(x)max(repmat(u0,1,params.nquad),w0(:,1)*x'+repmat(w0(:,2),1,params.nquad));
    %a. Calculate quadratures.
    quad_firm = jhandle(pts)*wgt;
    quad_worker = whandle(pts)*wgt;
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
    %To threshold pi's, identify zero and nonzero slope coefficients of job and work valfuns
    pi_h_firm = zeros(params.nshocks,1); pi_h_worker = zeros(params.nshocks,1);
    nonzero_j = (j0(:,1)~=0); zero_j = (j0(:,1)==0);
    nonzero_w = (w0(:,1)~=0); zero_w = (w0(:,1)==0);
    %For nonzero coefficients a.
    pi_h_firm(nonzero_j,1) = (v0(nonzero_j)-j0(nonzero_j,2))./j0(nonzero_j,1);
    pi_h_worker(nonzero_w,1) = (u0(nonzero_w)-w0(nonzero_w,2))./w0(nonzero_w,1);
    % Part 2. Given threshold values, compute next period J and W at
    %pi = 1. By constuction, we know that J(pi= pi_hat) = V.
    J1 = zprod.*params.yg - feval(wagefun,1) + params.beta*(1-params.lambda)*transmat*sum(j0,2);
    W1 = feval(wagefun,1) - params.adisutil + params.beta*transmat*(params.lambda*u0+(1-params.lambda)*sum(w0,2));
    Jpihf = v0;
    Wpihw = u0;
    %Second, compute new coefficients of linear functions.
    job(:,1) = (J1-Jpihf)./(1-pi_h_firm);
    job(:,2) = Jpihf - job(:,1).*pi_h_firm;
    work(:,1) = (W1-Wpihw)./(1-pi_h_worker);
    work(:,2) = Wpihw - work(:,1).*pi_h_worker;
    
    %Compare updated coefficients with old ones.
    delta_j = abs(job-j0); delta_v = abs(vac-v0);
    delta_w = abs(work-w0); delta_u = abs(unemp-u0);
    
    %Compute convergence criterion.
    delta = sum(sum(delta_j+delta_w)) + sum(delta_v+delta_u);
    
    %Update the guess.
    v0 = vac; u0 = unemp;
    j0 = job; w0 = work;  
end
toc
counter
end