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
h = 1/128;  %minimun gird spacing
num = fix(2 * RInitShrink / h);
r = linspace(0, 2 * RInitShrink, num);

epi = 0.005;
gamma = 10;

dt = 0.0001;

%variables

%Volumetric Compression Growth
SVcGrow = (1 - 2 * v) * Fz / (E * pi * RInitGrow^2);
cVcGrow = SVcGrow / k;
VlinVcGrow = alphaVcGrow * cVcGrow - beta;

%Strain Energy Density Growth
SSedGrow = Fz^2 / (2 * E * pi^2 * RInitGrow^4);
cSedGrow = SSedGrow / k;
VlinSedGrow = alphaSedGrow * cSedGrow - beta;

%Volumetric Compression Shrinkage
SVcShrink = (1 - 2 * v) * Fz / (E * pi * RInitShrink^2);
cVcShrink = SVcShrink / k;
VlinVcShrink = alphaVcShrink * cVcShrink - beta;

%Strain Energy Density Shrinkage
SSedShrink = Fz^2 / (2 * E * pi^2 * RInitShrink^4);
cSedShrink = SSedShrink / k;
VlinSedShrink = alphaSedShrink * cSedShrink - beta;


%% 1d Growth - Volumetric Compression
%%% working


phi = zeros(20000,num);
phi(:,1) = 1;
phi(:,num) = 0;

phi(1,:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

plot(r,phi(1,:));
hold on;

for j = 2:20000
    for i = 2:num-1
        phi(j,i) = phi(j-1,i) + dt * (-VlinVcGrow * (phi(j-1,i) - phi(j-1,i-1))/h + ...(-phi(j-1,i) * (1 - phi(j-1,i)) / epi) + ...
            gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i)) + ...
            gamma * epi^2 * (phi(j-1,i-1) - 2 * phi(j-1,i) + phi(j-1,i+1))/h^2);...(phi(j-1,i) - phi(j-1,i-1))^2 / h^2);    
    end
    if mod(j,200) == 0
        plot(r,phi(j,:));
    end
    
end

%% Working with funtions in SOLIDIFICATION

% phi = zeros(20000,num);
% phi(:,1) = 1;
% phi(:,num) = 0;
% 
% phi(1,:) = (1 -tanh((r - RInitGrow)/2/epi))/2;
% plot(r,phi(1,:));
% hold on;
% 
% 
% for j = 1 : 20000
%     for i = 2 : num-1
%         phi(j+1,i) = phi(j,i) + dt * (-VlinSedGrow * (-phi(j,i) * (1 - phi(j,i)) / epi) + ...
%             gamma * epi^2 * ((phi(j,i-1) - 2 * phi(j,i) + phi(j,i+1)) / h^2 + ...
%             phi(j,i) * (1 - phi(j,i)) * (1 - 2 * phi(j,i)) / epi^2 ));
%     end
%     if mod(j,200) == 0
%         plot(r,phi(j,:));
%     end
% end
% plot(r,phi(200,:));
% 

%% Testing with funtions in SHARP INTERFACE
% Working but wrong direction

% phi = zeros(200,64);
% phi(:,1) = 1;
% phi(:,64) = -1;
% 
% phi(1,:) = -tanh(r/sqrt(2)/epi);
% plot(r,phi(1,:));
% hold on;
% 
% for j = 1 : 199
%     for i = 2 : 63
%         phi(j+1,i) = phi(j,i) + dt * (-VlinVcGrow * (1 - phi(j,i)^2) / sqrt(2) / epi + ...
%             gamma * epi^2 * ((phi(j,i) - phi(j,i-1))^2 / h^2 + ...
%             phi(j,i) * (1 - phi(j,i)^2) / epi^2 ));
%     end
% end
% plot(r,phi(200,:));