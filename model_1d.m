%parameters
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

epi = 0.005;
gamma = 10;

dt = 0.01;

%variables

%Volumetric Compression Growth
SVcGrow = (1 - 2 * v) * Fz / (E * pi * RInitGrow^2);
cVcGrow = SVcGrow / k;
VlinVcGrow = alphaVcGrow * cVcGrow - beta;

%Strain Energy Density Growth
SSedGrow = Fz^2 / (2 * E * pi^2 * RInitGrow^4);
cSedGrwo = SSedGrow / k;
VlinSedGrow = alphaSedGrow * cSedGrow - beta;

%Volumetric Compression Shrinkage
SVcShrink = (1 - 2 * v) * Fz / (E * pi * RInitShrink^2);
cVcShrink = SVcShrink / k;
VlinVcShrink = alphaVcShrink * cVcShrink - beta;

%Strain Energy Density Shrinkage
SSedShrink = Fz^2 / (2 * E * pi^2 * RInitShrink^4);
cSedShrink = SSedShrink / k;
VlinSedShrink = alphaSedShrink * cSedShrink - beta;


%begin
r = linspace(-RInitGrow,RInitGrow,2*RInitGrow/h);

phi = zeros(200,64);
phi(:,1) = 1;
phi(:,64) = 0;

phiInit = 0.5 * (1 - tanh(r/(sqrt(8) * epi)));

phi(1,:) = phiInit;

plot(r,phi(1,:));
hold on;
for j = 2:200
    for i = 2:63
        phi(j,i) = phi(j-1,i) + dt * (- VlinVcGrow * abs((phi(j-1,i+1)-phi(j-1,i)) / h) + ...
            gamma * (-phi(j-1,i)^3 + 1.5 * phi(j-1,i)^2 - 0.5 * phi(j-1,i)) + ...
            gamma * epi^2 * (phi(j-1,i-1) - 2 * phi(j-1,i) + phi(j-1,i+1)) / h^2);
    end
%     if mod(j,5) == 0
%         plot(r,phi(j,:));
%     end
    
end



% phi2 = -tanh(r/sqrt(2)/epi);
% plot(r,phi2);
% hold on;
% 
% for i = 2 : 63
%     phi1(i) = phiInit(i) + dt * (gamma * epi^2 * (phiInit(i-1) -2 * phiInit(i) + phiInit(i+1)) / h^2 + ...
%         gamma * phiInit(i) * (1 - phiInit(i)^2) - VlinSedGrow * abs(phiInit(i+1) - phiInit(i)) / h);
% end
% plot(r,phi1);