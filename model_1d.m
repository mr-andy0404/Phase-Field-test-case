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


%variables

%Volumetric Compression Growth
SVcGrow = (1 - 2 * v) * Fz / (E * pi * RInitGrow^2);
cVcGrow = SVcGrow / k;
VlinvsGrow = alphaVcGrow * SVcGrow / k - beta;

%Strain Energy Density Growth
SSedGrow = Fz^2 / (2 * E * pi^2 * RInitGrow^4);
cSedGrwo = SSedGrow / k;
VlinSedGrow = alphaSedGrow * SSedGrow / k - beta;

%Volumetric Compression Shrinkage
SVcShrink = (1 - 2 * v) * Fz / (E * pi * RInitShrink^2);
cVcShrink = SVcShrink / k;
VlinVcShrink = alphaVcShrink * SVcShrink / k - beta;

%Strain Energy Density Shrinkage
SSedShrink = Fz^2 / (2 * E * pi^2 * RInitShrink^4);
cSedShrink = SSedShrink / k;
VlinSedShrink = alphaSedShrink * SSedShrink / k - beta;


%begin
r = linspace(-RInitGrow,RInitGrow,2*RInitGrow/h);

phiInit = 0.5 * (1 - tanh(r/(sqrt(8) * epi)));

nabla = gradient(phiInit);

lap = gradient(nabla);

phi1 = zeros(1,64);
phi1(1) = 1;
phi1(64) = 0;

for i = 2:63
    phi1(i) = 0.01 * (phiInit(i)