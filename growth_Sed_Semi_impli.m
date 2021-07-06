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
h = 1/1280;  %minimun gird spacing
num = fix(2 * RInitShrink / h);
r = linspace(0, 2 * RInitShrink, num);

epi = 0.005;
gamma = 10;

dt = 0.001;

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

%% 1d Growth - Strain Energy Density
% working
% problem with |nablaphi|

phi = zeros(2/dt,num);
phi(:,1) = 1;
phi(:,num) = 0;

phi(1,:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

plot(r,phi(1,:));
hold on;
title('\phi changing with time using Sed');
xlabel('r');
ylabel('\phi');
ylim([0 1]);

R = zeros(1,2/dt);
R(1) = RInitGrow;
A = zeros(num-2,num-2);
B = zeros(num-2,1);

for j = 2:2/dt
    
    %%updating velocity
    SSedGrow = Fz^2 / (2 * E * pi^2 * R(j-1)^4);
    cSedGrow = SSedGrow/k;
    VlinSedGrow = alphaSedGrow * cSedGrow - beta;
   
    for i = 2:num-1

        nablaphi = -2 * phi(j-1,i) * (1 - phi(j-1,i)) / sqrt(8) / epi;
        lapphi = -phi(j-1,i) * (1 - phi(j-1,i)) * (1 - 2 * phi(j-1,i)) / 2 / epi^2;
        
        A(i-1,i-1) = 1 + 2 * dt * gamma * epi^2 / h^2;
        if i < num-1
            A(i-1,i) = - dt * gamma * epi^2 / h^2;
        end
        if i > 2
            A(i-1,i-2) = - dt * gamma * epi^2 / h^2;
        end
        
        if i == 2
            B(i-1) = dt * gamma * epi^2 / h^2 + phi(j-1,i) + phi(j-1,i) + ...
                dt * (-VlinSedGrow * nablaphi + ...
                gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i))); 
        else
            B(i-1) = phi(j-1,i) + dt * (-VlinSedGrow * nablaphi + ...
            gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i)));
        end
         
    end
    
    phi(j,2:num-1) = A\B;
    
     R(j) = 0.5 * (max(r(phi(j,:) >= 0.5)) + min(r(phi(j,:) < 0.5)));
    
    if mod(j,dt*1e4) == 0
        plot(r,phi(j,:));
    end
    
end

figure(2);
t = linspace(1,2/dt,2/dt)* dt;
plot(t, pi * R.^2);
title('Volume varying with time with Sed');
xlabel('time/s');
ylabel('Volume');
ylim([0.15,0.65]);
grid on;

