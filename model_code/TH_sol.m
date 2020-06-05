% the last ~40 lines of gc_molND_nonuniform_s.m
H0 = 0.1;
E = 5;
E1 = 1;
E2 = 1;
E3 = 2;
E4 = 3;
qbarUA = 1;  % initial volume for unit width channel (see integral of 
                    % initial condition for the case theta = 0 which gives
                    % qbarUA = 2 * H0 / 3
B_TH = beta(E1/E3,1 + E2/E4);
eta_N = (E3 * ( (E*E3)/(E2*E4) )^(E2/E4) /B_TH  )^(E4/E);
num_eta = n;
eta = s*eta_N;
phi = ( (E2*E4)/(E*E3) * (eta_N^E3 - eta.^E3) ).^(1/E4);

t_TH_f = t(end);
xN_TH_f = eta_N* ( (1/3)*qbarUA^E4 * t_TH_f).^(1/E);
h_TH_f = (qbarUA^2 * 3./t_TH_f).^(1/E) * phi;