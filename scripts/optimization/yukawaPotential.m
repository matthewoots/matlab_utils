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

data = linspace(0,10,100);
d = data;
% Increasing alpha shrinks the area
alpha = 0.3;
A = 1;
boundary_val = -20;
% should we add a (-) to the denominator :
% the whole curve is flipped
for i=1:numel(d)
    J(i) = A * (exp(-alpha * d(i))) / -d(i);
    if isinf(J(i))
        J(i) = boundary_val;
    end
end
figure(1)
hold on
plot(d,J,'x');
grid on;