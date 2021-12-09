%% Yukawa Potential 
% Taken from the quantum field theory, 
% 	describing the potential of subatomic particles

% Hence particles right in the center of each other ranges to inf

% A = Scaling factor 
% alpha = decay rate
% d = cost factor
% J = A * (e ^ (-alpha * d) / d);
clc
clear all
close all 

% data can only be positive

data = linspace(0.0,10,100);
d = data;
% Increasing alpha shrinks the area
alpha = 0.1;
A = 1;
boundary_val = 0.05;
% should we add a (-) to the denominator :
% the whole curve is flipped
offset = (exp(-alpha * boundary_val)) / boundary_val;
for i=1:numel(d)
    J0(i) = (exp(-alpha * d(i))) / -d(i) + offset;
    J1(i) = (exp(-(alpha*10) * d(i))) / -d(i) + offset;
    J2(i) = (exp(-(alpha*100) * d(i))) / -d(i) + offset;
    if isinf(J0(i))
        J0(i) = 0;
    end
    if isinf(J1(i))
        J1(i) = 0;
    end
    if isinf(J2(i))
        J2(i) = 0;
    end
end
figure(1)
hold on
plot(d,J0,'x','DisplayName',string(alpha));
plot(d,J1,'x','DisplayName',string(alpha*10));
plot(d,J2,'x','DisplayName',string(alpha*100));
legend
xlabel('distance')
ylabel('cost');
grid on;