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


%% 1d Growth - Strain Energy Density
% working

phi = zeros(2/dt,num);
phi(:,1) = 1;
phi(:,num) = 0;
gradphi = zeros(1,num);

phi(1,:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

plot(r,phi(1,:));
hold on;
title('\phi changing with time using Vc');
xlabel('r');
ylabel('\phi');
% ylim([0 1]);

R = zeros(1,2/dt);
R(1) = RInitGrow;
A = zeros(num-2,num-2);
B = zeros(num-2,1);

for j = 2:2/dt
    
    %%updating velocity
    SVcGrow = (1 - 2 * v) * Fz / (E * pi * R(j-1)^2);
    cVcGrow = SVcGrow / k;
    VlinVcGrow = alphaVcGrow * cVcGrow - beta;
    gradphi = -2 * phi(j-1,:) .* (1 - phi(j-1,:)) / sqrt(8) / epi;

    for i = 2:num-1
        
%         nablaphi = -2 * phi(j-1,i) * (1 - phi(j-1,i)) / sqrt(8) / epi;
        
        A(i-1,i-1) = 1 + 2 * dt * gamma * epi^2 / h^2;
        if i < num-1
            A(i-1,i) = - dt * gamma * epi^2 / h^2;
        end
        if i > 2
            A(i-1,i-2) = - dt * gamma * epi^2 / h^2;
        end
        
        if i == 2
            B(i-1) = dt * gamma * epi^2 / h^2 + phi(j-1,i) + ...
                dt * (-VlinVcGrow * gradphi(i) + ...
                gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i))); 
        else
            B(i-1) = phi(j-1,i) + dt * (-VlinVcGrow * gradphi(i) + ...
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
title('Volume varying with time with Vc');
xlabel('time/s');
ylabel('Volume');
ylim([0.15,0.65]);
grid on;

