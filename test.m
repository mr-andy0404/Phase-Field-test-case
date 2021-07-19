d = 1;
S = 1.7;
k = 1/2;

RInitGrow = sqrt(0.2/pi);
RInitShrink = sqrt(0.6/pi);

h = 1/1280;  %minimun gird spacing
num = fix(2 * RInitShrink / h);
r = linspace(0, 2 * RInitShrink, num);
epi = 0.005;
gamma = 10;

dt = 0.001;

phi = zeros(1,num);
phi(1) = 1;
phi(num) = 0;
d2 = 1e-3;
k2 = 1e3;

phi(:) = 0.5 * (1 - tanh((r-RInitGrow)/(sqrt(8) * epi)));

CA = zeros(num-2,num-2);
CB = zeros(num-2,1);

% for i = 1:num-2 % with c(0) = 3.4 and c(r) = 1.1333
%     CA(i,i) = -2/h^2;
%     if i > 1
%         CA(i,i-1) = 1/h^2;
%     end
%     if i < num-2
%         CA(i,i+1) = 1/h^2;
%     end
%     
%     if i == 1
%         CB(i) = -S/d - 3.4/h^2;
%     elseif i == num-2
%         CB(i) = -S/d - 1.1333/h^2;
%     else
%         CB(i) = -S/d;
%     end
% end

% for i = 1:num-2
%     CA(i,i) = -2 * (d * phi(i+1) + d2 * (1 - phi(i+1))) /h^2;
%     if i > 1
%         CA(i,i-1) = (d * phi(i+1) + d2 * (1 - phi(i+1)))/h^2;
%     end
%     if i < num-2
%         CA(i,i+1) = (d * phi(i+1) + d2 * (1 - phi(i+1)))/h^2;
%     end
%     
%     if i == 1
%         CB(i) = -S * phi(i+1) - 3.4 * (d * phi(i+1) + d2 * (1 - phi(i+1)))/h^2;
%     elseif i == num-2
%         CB(i) = -S * phi(i+1) - 1.1333 * (d * phi(i+1) + d2 * (1 - phi(i+1)))/h^2;
%     else
%         CB(i) = -S * phi(i+1);
%     end
% end

for i = 1:num
    CA(i,i) = 2 * (d * phi(i) + d2 * (1 - phi(i))) /h^2;
    if i > 1
        CA(i,i-1) = -(d * phi(i) + d2 * (1 - phi(i)))/h^2;
    end
    if i < num
        CA(i,i+1) = -(d * phi(i) + d2 * (1 - phi(i)))/h^2;
    end
    
        CB(i) = S * phi(i);

end


c = CA\CB;

plot(r,c);




