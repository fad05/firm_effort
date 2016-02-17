      %Given the value functions, find vector of threshold values of pi.
    %To do it, identify zero and nonzero slopes of job and work valfuns
    nonzero_j = (j0(:,1)~=0); zero_j = (j0(:,1)==0);
    nonzero_w = (w0(:,1)~=0); zero_w = (w0(:,1)==0);
    %For nonzero coefficients a.
    pi_h_firm(nonzero_j,1) = (v0(nonzero_j)-j0(nonzero_j,2))./j0(nonzero_j,1);
    pi_h_worker(nonzero_w,1) = (u0(nonzero_w)-w0(nonzero_w,2))./w0(nonzero_w,1);
    pi_h_firm(pi_h_firm>1) = 1;
    pi_h_firm(pi_h_firm<0) = 0;
    pi_h_worker(pi_h_worker>1) = 1;
    pi_h_worker(pi_h_worker<0) = 0;

%Second, for zero elements.
    pi_h_firm(zero_j & j0(zero_j,2)>=v0(zero_j),1) = 0;
    pi_h_firm(zero_j & j0(zero_j,2)<v0(zero_j),1) = 1;
    pi_h_worker(zero_w & w0(zero_w,2)>=u0(zero_w),1) = 0;
    pi_h_worker(zero_w & w0(zero_w,2)<u0(zero_w),1) = 1;