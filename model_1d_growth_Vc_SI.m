clear;

%% parameters
E = 1;      %Young's modulus
v = 0.33;   %Poisson ratio

RInitGrow = sqrt(0.2/pi);
RInitShrink = sqrt(0.6/pi);

alphaSedGrow = 0.36;
alphaVcGrow = 0.9;
alphaSedShrink = 0.04;
alphaVcShrink = 0.3;

beta = 1;

Fz = 1;
d = 1;
k = 1/2;
h = 1/10000;  %minimun gird spacing
num = fix(2 * RInitShrink / h);
r = linspace(-RInitShrink, RInitShrink, num);

epi = 0.005;
gamma = 10;

dt = 0.0001;

%% variables

% %Volumetric Compression Growth
% SVcGrow = (1 - 2 * v) * Fz / (E * pi * RInitGrow^2);
% cVcGrow = SVcGrow / k;
% VlinVcGrow = alphaVcGrow * cVcGrow - beta;
% 
% %Strain Energy Density Growth
% SSedGrow = Fz^2 / (2 * E * pi^2 * RInitGrow^4);
% cSedGrow = SSedGrow / k;
% VlinSedGrow = alphaSedGrow * cSedGrow - beta;
% 
% %Volumetric Compression Shrinkage
% SVcShrink = (1 - 2 * v) * Fz / (E * pi * RInitShrink^2);
% cVcShrink = SVcShrink / k;
% VlinVcShrink = alphaVcShrink * cVcShrink - beta;
% 
% %Strain Energy Density Shrinkage
% SSedShrink = Fz^2 / (2 * E * pi^2 * RInitShrink^4);
% cSedShrink = SSedShrink / k;
% VlinSedShrink = alphaSedShrink * cSedShrink - beta;

%% 1d Growth - Volumetric Compression - Sharp Interface

phi = zeros(20000,num);
phi(:,1) = 1;
phi(:,num) = -1;

phi(1,:) = - tanh((r) / sqrt(2) /epi); %Origin from -1 to 1, or change to 0 to 1?

plot(r,phi(1,:));
hold on;
title('\phi changing with time using Vc');
xlabel('r');
ylabel('\phi');
%ylim([0 1]);

R = zeros(1,20000);
R(1) = RInitGrow;

for j = 2:20000
    
    %%updating velocity
    SVcGrow = (1 - 2 * v) * Fz / (E * pi * R(j-1)^2);
    cVcGrow = SVcGrow / k;
    VlinVcGrow = alphaVcGrow * cVcGrow - beta;
% VlinVcGrow = 1;
    
    for i = 2:num-1
        
        nablaphi = (1 - phi(j-1,i)^2) / sqrt(2) / epi;
        lapphi = phi(j-1,i) * (1 - phi(j-1,i)^2) * (1 - 2 * phi(j-1,i)) / 2 / epi^2;
        
        phi(j,i) = phi(j-1,i) + dt * (-VlinVcGrow * nablaphi + ...
            gamma * phi(j-1,i) * (1 - phi(j-1,i)^2) + ...
            gamma * epi^2 * nablaphi^2);  
        
    end
    
    R(j) = 0.5 * (max(r(phi(j,:) >= 0.5)) + min(r(phi(j,:) < 0.5)));
    
    if mod(j,200) == 0
        plot(r,phi(j,:));
    end
    
end

figure(2);
t = linspace(1,20000,20000)* dt;
plot(t, pi * R.^2);
title('Volume varying with time using Vc');
xlabel('time/s');
ylabel('Volume');
%ylim([0.15,0.65]);
grid on;

