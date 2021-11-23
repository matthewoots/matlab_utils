%% Fixed Wing
% A fixed wing aircraft dynamic can be expressed by a point-mass mode
% [1] H. Alturbeh and J. F. Whidborne. Real-time obstacle collision avoidance for fixed wing aircraft using b-splines. In 2014 UKACC International Conference on Control (CONTROL), pages 115–120, 2014.

% Simplified B-spline calculation using Beinstein basis formula 
%   implemented from https://github.com/felaube/b-spline-trajectories
%   to B-Spline in matrix representation

% gamma = Incline angle
% phi = Heading angle

clc
clear all 
close all

funct_path = '../../functions';
addpath(strcat(funct_path,'/common'));
addpath(strcat(funct_path,'/bspline'));

% Gravitational acceleration [m/s²]
g = 9.81;
% Air density (speed is considered to be AES) [kg/m³]
rho = 1.225;

%% Aircraft Parameters
% Oswald factor [-]
e = 0.8;

% Wing area [m²]
S = 1.554;

% Span [m]
b = 3.9422;

% Maximum thrust (all operating engines) [N]
T_max = 130;

% Aircraft mass [kg]
mass = 25;

% Aspect Ratio [-]
AR = b^2/S;

% k value [-]
k = 0.0334;

% Cd0 value [-]
C_D_0 = 0.0715;

% Bank Angle min and max
gamma_max = deg2rad(20);
gamma_min = deg2rad(-20);

% Drag minimum
V_min_drag = sqrt(mass * g * 2 / (rho*S))*((k/C_D_0)^(1/4));

% Minumum trust
T_min = 1/2 * rho * S * 2 * C_D_0 .* V_min_drag .^2;

%% Use bspline to get uav states
order = 4;
k = order + 1;

M = getM(order);
u = zeros(1,k); du = zeros(1,k); ddu = zeros(1,k); dddu = zeros(1,k);
p = zeros(k,1);

% p = degree
% n = cp - 1
% knots = m + 1
% basics is m = n + p + 1
% Open B-splines are [knot(p), knot(m-p)]
% * Clamping only happens for the knots not for the control points!
% * Control points are not really attached to time

% User inputs
timespan = [0,20];
pos = []; vel = []; acc = []; jrk = []; tt = [];
kn_seg = 10; % Division of 1 span or 1 segment
for i=1:3
    r = rand(1,4);
    for j = 1:order
        r = [r(1) r r(end)];
    end
    numrange = 12; % range for evaluation
    P(i,:) = numrange * r; % multiply the random [0-1] by the number range
    n = numel(r) - 1; % n = controlpoints - 1

    % How much to segment a knotspan (internal evaluation in a spline segment)
    buffer = zeros(1,order + 1);
    knots_total = n + order + 1 + 1;
    % Range of index to evaluate according
    range = [order:numel(P(i,:))]; 
    % range = [order:knots_total-order]; 
    % Not <---order---><----Accepted Range numel(P)-----><----order+1---->
    % But <---order---><----------Accepted Range numel(P)---------------->

    % I may have made the code abit different by not using knots = m + 1 which
    % means using numel(P) rather than knots_total to segment the time horizon
    % Solution for matching time span
    dt = (timespan(2) - timespan(1)) / (numel(range)-1); % span of 1 knot
    t = linspace(timespan(1), timespan(2), (numel(range))); % knots
    % End of solution for matching time span

    [p0,v0,a0,t0] = getBSpline(order, timespan, P(i,:), kn_seg, false);
    pos(i,:) = p0; 
    vel(i,:) = v0; 
    t = t0;
end

figure(1)
hold on
for i = 1:numel(range)-1
    % t(i) should be empty and the control point matching with the time
    % should be +1, that is why t(i+1)
    plot3(P(1,range(i)),P(2,range(i)),P(3,range(i)),'o');
end
plot3(pos(1,:),pos(2,:),pos(3,:),'x');
grid on

figure(2)
subplot(2,1,1)
hold on
plot(t,pos(1,:))
plot(t,pos(2,:))
plot(t,pos(3,:))
title("Position")
hold off
grid on

subplot(2,1,2)
hold on
plot(t,vel(1,:))
plot(t,vel(2,:))
plot(t,vel(3,:))
title("Velocity")
hold off
grid on

%% Fixed wing optimization
% Calculate flight path angle
% gamma_a = asin(w_a./V_a);

% heading angle and their derivatives