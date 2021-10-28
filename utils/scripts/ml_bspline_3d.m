clc
clear all
close all

funct_path = '../functions/bspline_utils';
addpath(funct_path);

% Illustrates B-spline curve estimation without knowing parameter values.
% Copyright 2010 Levente Hunyadi
% spline order
k = 4;
% knot sequence
t_0 = [0 0 0 0 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1 1 1];

% control points (unknown)
D_0 = [ 0.1993 0.4965 0.6671 0.7085 0.6809 0.2938 0.1071 0.3929 0.5933 0.8099 0.8998 0.8906 ...
      ; 0.8377 0.8436 0.7617 0.6126 0.212 0.1067 0.3962 0.5249 0.5015 0.3991 0.6477 0.8553 ...
      ; 0.1993 0.4965 0.6671 0.6126 0.212 0.1067 0.3962 0.5249 0.5015 0.3991 0.6477 0.8553 ];
tic
% points on B-spline curve
p = bspline_deboor(k,t_0,D_0);
pt = linspace(t_0(end),t_0(1),length(p));
fprintf("p optimization time %f\n",toc);

tic
% velocity on B-spline curve
[t_1,D_1] = bspline_deriv(k-1, t_0, D_0);
v = bspline_deboor(k-1,t_1,D_1);
fprintf("v optimization time %f\n",toc);
vt = linspace(t_1(end),t_1(1),length(v));

% plot control points and spline
figure;
subplot(2,2,1:2)
hold on;
plot3(p(1,:), p(2,:), p(3,:), 'x');
% plot control points
for i=1:length(D_0)
    plot3(D_0(1,i), D_0(2,i), D_0(3,i), 'o');
end
% plot line segments
for i=1:length(D_0)-1
    plot3(D_0(1,i:i+1), D_0(2,i:i+1), D_0(3,i:i+1), '-');
end
hold off;
legend('original curve');
title('Trajectory');
grid on

% Plot position
subplot(2,2,3)
hold on;
plot(pt, p(1,:)); plot(pt, p(2,:)); plot(pt, p(3,:));
legend('x','y','z');
title('Position graph');
grid on
hold off;

% Plot velocity
subplot(2,2,4)
hold on;
plot(vt, v(1,:)); plot(vt, v(2,:)); plot(vt, v(3,:));
legend('vx','vy','vz');
title('Velocity graph');
grid on
hold off;
