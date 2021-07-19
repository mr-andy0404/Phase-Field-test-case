clear;

%% parameters
E = 1;      %Young's modulus
v = 1/3;   %Poisson ratio

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
h = 1/10000;  %minimun gird spacing
num = fix(2 * RInitShrink / h);
r = linspace(0, 2 * RInitShrink, num);

epi = 0.005;
gamma = 10;

dt = 0.0001;


%% 1d Growth - Volumetric Compression

phi = zeros(20000,num);
phi(:,1) = 1;
phi(:,num) = 0;
gradphi = zeros(1,num);

phi(1,:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

plot(r,phi(1,:));
hold on;
title('\phi changing with time');
xlabel('r');
ylabel('\phi');

RVc = zeros(1,20000);
RVc(1) = RInitGrow;
AVc = zeros(num-2,num-2);
BVc = zeros(num-2,1);

for j = 2:20000
    
 %%updating velocity
    SVcGrow = (1 - 2 * v) * Fz / (E * pi * RVc(j-1)^2);
    cVcGrow = SVcGrow / k;
    VlinVcGrow = alphaVcGrow * cVcGrow - beta;
    gradphi = -2 * phi(j-1,:) .* (1 - phi(j-1,:)) / sqrt(8) / epi;

    for i = 2:num-1
        
        
        AVc(i-1,i-1) = 1 + 2 * dt * gamma * epi^2 / h^2;
        if i < num-1
            AVc(i-1,i) = - dt * gamma * epi^2 / h^2;
        end
        if i > 2
            AVc(i-1,i-2) = - dt * gamma * epi^2 / h^2;
        end
        
        if i == 2
            BVc(i-1) = dt * gamma * epi^2 / h^2 + phi(j-1,i) + ...
                dt * (-VlinVcGrow * gradphi(i) + ...
                gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i))); 
        else
            BVc(i-1) = phi(j-1,i) + dt * (-VlinVcGrow * gradphi(i) + ...
            gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i)));
        end
         
    end
    
    phi(j,2:num-1) = AVc\BVc;
    
     RVc(j) = 0.5 * (max(r(phi(j,:) >= 0.5)) + min(r(phi(j,:) < 0.5)));
    
    if mod(j,dt*1e4) == 0
        plot(r,phi(j,:));
    end
end

figure(2);
t = linspace(1,20000,20000)* dt;
plot(t, pi * RVc.^2);
title('Volume varying with time');
xlabel('time/s');
ylabel('Volume');
ylim([0.15,0.65]);
grid on;
hold on;


%% 1d Growth - Strain Energy Density
%%% working

phi = zeros(20000,num);
phi(:,1) = 1;
phi(:,num) = 0;

phi(1,:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

figure(3);
plot(r,phi(1,:));
hold on;
title('\phi changing with time');
xlabel('r');
ylabel('\phi');

RSed = zeros(1,20000);
RSed(1) = RInitGrow;

for j = 2:20000

        nablaphi = -2 * phi(j-1,i) * (1 - phi(j-1,i)) / sqrt(8) / epi;
        lapphi = -phi(j-1,i) * (1 - phi(j-1,i)) * (1 - 2 * phi(j-1,i)) / 2 / epi^2;
        
        A(i-1,i-1) = 1 + 2 * dt * gamma * epi^2 / h^2; % define the digonal element
        
        if i < num-1    % define the upper part
            A(i-1,i) = - dt * gamma * epi^2 / h^2;
        end
        if i > 2    % define the lower part
            A(i-1,i-2) = - dt * gamma * epi^2 / h^2;
        end
        
        
        if i == 2
            B(i-1) = dt * gamma * epi^2 / h^2 + phi(j-1,i) + ...
                dt * (-VlinSedGrow * nablaphi + ...
                gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i))); 
        else
            B(i-1) = phi(j-1,i) + dt * (-VlinSedGrow * nablaphi + ...
            gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i)));
        end
         
    
    phi(j,2:num-1) = A\B;
    
     R(j) = 0.5 * (max(r(phi(j,:) >= 0.5)) + min(r(phi(j,:) < 0.5)));
    
    if mod(j,dt*1e4) == 0
        plot(r,phi(j,:));
    end
 
end

figure(2);
t = linspace(1,20000,20000)* dt;
plot(t, pi * RSed.^2);
title('Volume varying with time');
xlabel('time/s');
ylabel('Volume');
ylim([0.15,0.65]);
legend('1d growth with Vc','1d growth with Sed');

