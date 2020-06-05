% close all;
% clear all;
% n=50;  %800; % >=4

% boundary interval
dr = 1/(n-1);
r = linspace(0,1,n);
% lambda = 5.5;
s = (1 - exp(-lambda*r))/(1 - exp(-lambda));

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ORIGINAL ATTEMPT AT AN INITIAL CONDITION
%
% % % H0 = 0.1;   % dimensionless curly H (see 8-14-19 notes, p.9 ... can relate to volume if we want)
% % % theta = 0.0;
% % % hIn = (H0 + tan(theta))*(1-s(1:n-1).^2) + ...
% % %       (s(1:n-1)-1)*tan(theta);
% % % yIn = [hIn 1]; %IC for h and IC for xN
% % % tInit  = 0;
% % % tFinal = 100000.0; %100000000.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SIMILARITY INITIAL CONDITION
%
theta = 0;
tInit = 1;
% tFinal = 10.0;
%
% match initial condition length with length of s(1:n-1)
%
E = 5;
E1 = 1;
E2 = 1;
E3 = 2;
E4 = 3;
qbarUA = 1; %(2*H0/3);  % initial volume for unit width channel (see integral of 
                    % initial condition for the case theta = 0 which gives
                    % qbarUA = 2 * H0 / 3
B_TH = beta(E1/E3,1 + E2/E4);
eta_N = (E3 * ( (E*E3)/(E2*E4) )^(E2/E4) /B_TH  )^(E4/E);
eta_INIT = s*eta_N;
phi_INIT = ( (E2*E4)/(E*E3) * (eta_N^E3 - eta_INIT(1:n-1).^E3) ).^(1/E4);
t_TH_INIT = tInit;
x_N_TH_INIT = eta_N * ( (1/3) * qbarUA^E4 * t_TH_INIT).^(1/E);
x_TH_INIT_plot = s*x_N_TH_INIT;
h_TH_INIT = (qbarUA^2 * 3./t_TH_INIT ).^(1/E) * phi_INIT;

yIn = [h_TH_INIT x_N_TH_INIT]; %IC for h and IC for xN

% plot the initial condition to help debug ... (DMA, 11-8-2019)
%%figure(88);plot(x_TH_INIT_plot,[h_TH_INIT 0]);hold on;
%%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% time domain, ode23 config
tspan = [tInit tFinal];
RelTolVal=10^(-3);
AbsTolVal=10^(-3);
options = odeset ('RelTol',RelTolVal,'AbsTol',AbsTolVal,'InitialStep',1e-3); %,'MaxStep',1e-8);

[t,y]  = ode23 (@gc_rhsND_nonuniform_s, tspan, yIn', options, n, r, dr, lambda, bta, theta, a);


%%%%
%%%% plotting
%%%%
% figure(1);plot(t,y(:,end))
% xlabel('time','FontSize',16);
% ylabel('xN','FontSize',16);
% %
% figure(4);loglog(t,y(:,end),'*');hold on;
% xlabel('time','FontSize',16);
% ylabel('xN','FontSize',16);
% %
% % estimate for power law xN = A t^B
% %
% B_est = log(y(end,end)/y(end-1,end))/log(t(end)/t(end-1));
% Ashift = 1;
% xNpowerlaw = Ashift * t.^B_est;
% loglog(t,xNpowerlaw,'r--');
% %
% %
% figure(2);
% plot(s,[y(1,1:n-1) 0]);hold on;
% %plot(s,[y(100,1:n-1) 0]);
% plot(s,[y(end,1:n-1) 0]);
% xlabel('s','FontSize',16);
% ylabel('h','FontSize',16);
% %
% 
% 
% figure(3);
% plot(s*y(1,end),[y(1,1:n-1) 0],'r--','LineWidth',3);hold on;   % initial condition
% ncurves = 6;
% for jj=1:ncurves;
%     curve_time_index = floor(jj/(ncurves+1) * length(t));
%     plot(s*y(curve_time_index,end),[y(curve_time_index,1:n-1) 0],'b-.','LineWidth',1);hold on;
% end
% plot(s*y(end,end),[y(end,1:n-1) 0],'k-','LineWidth',3);
% xlabel('x','FontSize',16);
% ylabel('h','FontSize',16);
% 
% 
% %%%%% Takagi/Huppert 2008 similarity solution
% %%%%% constant volume (a=0)
% %%%%% horizontal flow (H=1)
% %
% 
% %
% num_eta = n;
% eta = linspace(0,eta_N,num_eta);
% phi = ( (E2*E4)/(E*E3) * (eta_N^E3 - eta.^E3) ).^(1/E4);
% %
% % use the same choices for time as in Figure 3 to plot this similarity solution
% %
% figure(3)
% ncurves = 6;
% for jj=1:ncurves;
%     curve_time_index = floor(jj/(ncurves+1) * length(t));
%     t_TH_plot = t(curve_time_index);
%     x_N_TH = eta_N * ( (1/3) * qbarUA^E4 * t_TH_plot).^(1/E);
%     x_TH = linspace(0,x_N_TH,num_eta);
%     h_TH = (qbarUA^2 * 3./t_TH_plot ).^(1/E) * phi;
% %%    plot(eta,h_TH,'c--','LineWidth',3);hold on;
%     plot(x_TH,h_TH,'c--','LineWidth',3);hold on;
% end
% t_TH_f = t(end);
% xN_TH_f = eta_N* ( (1/3)*qbarUA^E4 * t_TH_f).^(1/E);
% x_TH_f = linspace(0,xN_TH_f,num_eta);
% h_TH_f = (qbarUA^2 * 3./t_TH_f).^(1/E) * phi;
% plot(x_TH_f, h_TH_f, 'm-', 'LineWidth',3);
% 
% 
% figure(101)
% ncurves = 6;
% for jj=1:ncurves;
%     curve_time_index = floor(jj/(ncurves+1) * length(t));
%     t_TH_plot = t(curve_time_index);
%     x_N_TH = eta_N * ( (1/3) * qbarUA^E4 * t_TH_plot).^(1/E);
%     loglog(t_TH_plot,x_N_TH,'c*');hold on;
% end