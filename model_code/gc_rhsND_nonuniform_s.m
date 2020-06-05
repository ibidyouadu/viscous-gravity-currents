function rhs = gc_rhsND_nonuniform_s(t, yIn, n, r, dr, lambda, bta, theta, a);
%
% define constants
%a=0;
q=1; % make this equal to qbarUA!
%W=1;
%rho=1;
%nu=1;
%g=1;
%%% theta=0.1;  %pi/6; % defined in gc_mol.m
b0=bta;
b1=bta;
b2=bta;   % dimensionless betas (CURRENTLY USING FOR PARAMETER TESTING)

% y = [h_0, h(s=0), ... h(s=1-ds), h(s=1), xN]
y = [0 yIn(1:n-1)' 0 yIn(n)]; % h0=y(1) will be set shortly; h_n = y(n+1) = 0 (BC)
h = y(1:n+1); % h0, h1, ..., hn

% construct s from r; compute drds
% using g2 from Dalwadi et. al
%   r = g2(s; lambda) so r->g2, s->eta
s = (1 - exp(-lambda*r))/(1 - exp(-lambda));
drds = (1 - exp(-lambda))./(lambda * (1 - s*(1 - exp(-lambda))));

% define terms for the rhs equations, 
f = 1/3*h.^3 + b0*h.^2 + (b1+b2)*h;
if a==0
    y(1) = y(3);
else
    y(1) = y(3) + 2*dr*y(end) * (a*q*t^(a-1)/(cos(theta)*f(2)) - tan(theta))/drds(1);
end

h_j = y(2:n); % h1, h2, ..., h_(n-1)
h_jp = y(3:n+1); % h2, h3, ..., hn
h_jm = y(1:n-1); % h0, h1, ..., h_(n-2)

f(1) = 1/3*y(1).^3 + b0*y(1).^2 + (b1+b2)*y(1);


dhdr = (h_jp - h_jm)/(2*dr); % defined on h_1,...,h_(n-1)
dhdr_xN = (3*y(n+1) - 4*y(n) + y(n-1))/(2*dr);  % defined on h_n; for dxNdt eq

fp = (f(3:n+1) + f(2:n))/2; % don't need fp(h_0)
fm = (f(2:n) + f(1:n-1))/2; % don't need fm(h_n)

drds_p = (drds(2:n) + drds(1:n-1))/2;
drds_m = (drds(1:n-1) + [drds(1) drds(1:n-2)])/2;

Gp = fp .* (drds_p/y(end).*(h_jp - h_j)/dr - tan(theta)); % dont need Gp(h_0)
Gm = fm .* (drds_m/y(end).*(h_j - h_jm)/dr - tan(theta)); % dont need Gm(h_n)
dGdr = (Gp - Gm)/dr;

% define rhs
% this is for all n-1 dhdt and one dxNdt eq'ns
% indexing changes here; idx=1:n-1 is for h_1 to h_(n-1), idx=n is for xN
DY(n) =  - cos(theta)*(b1 + b2)*(drds(end)*dhdr_xN/y(end) - tan(theta));
DY(1:n-1) = s(1:n-1)/y(end)*DY(n).*dhdr.*drds(1:n-1) + cos(theta)*(dGdr/y(end).*drds(1:n-1));
rhs = DY';