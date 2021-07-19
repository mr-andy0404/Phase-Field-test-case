
clear;

%% parameters
E = 1;      %Young's modulus
v = 0.33;   %Poisson ratio

RInitGrow = sqrt(0.2/pi);
RInitShrink = sqrt(0.6/pi);

alphaSedGrow = 0.36;
alphaVcGrow = 0.8823;
alphaSedShrink = 0.04;
alphaVcShrink = 0.3;

beta = 1;

Fz = 1;
d = 1;
k = 1/2;
h = 1/1280;  %minimun gird spacing
num = fix(2 * RInitShrink / h);
r = linspace(0, 2 * RInitShrink, num);

epi = 0.005;
gamma = 10;

dt = 0.001;

phi = zeros(2/dt,num);
phi(:,1) = 1;
phi(:,num) = 0;

phi(1,:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

plot(r,phi(1,:));
hold on;
title('\phi changing with time using Vc');
xlabel('r');
ylabel('\phi');
% % ylim([0 1]);

R = zeros(1,2/dt);
R(1) = RInitGrow;
A = zeros(num-2,num-2);
B = zeros(num-2,1);
CA = zeros(num,num);
CB = zeros(num,1);
Ke = 1/h * [1,-1;-1,1];
Me = h/6 * [2,1;1,2];
Fe = h/2 * [1;1];

k2 = 1e3;
d2 = 1e-3;

    SVcGrow = (1 - 2 * v) * Fz / (E * pi * RInitGrow^2);
    cVcGrow = SVcGrow / k;
    VlinVcGrow = alphaVcGrow * cVcGrow - beta;
    gradphi = -2 * phi(1,:) .* (1 - phi(1,:)) / sqrt(8) / epi;

for ie = 1:num-1
    
    phim(ie) = (phi(1,ie)+phi(1,ie+1))/2;
    
    ke(ie) = k * phim(ie) + k2 * (1 - phim(ie));
    de(ie) = d * phim(ie) + d2 * (1 - phim(ie));
    
    CA(ie:ie+1,ie:ie+1) = CA(ie:ie+1,ie:ie+1) + ke(ie) * Me + de(ie) * Ke;
    
    CB(ie:ie+1) = CB(ie:ie+1) + SVcGrow * Fe * phim(ie);
    
    
end

    CVcGrow = CA\CB;
    VVcGrow = CVcGrow * alphaVcGrow - beta;
    plot(r,CVcGrow,r,VVcGrow);
    legend("phi","concentration","velocity");

