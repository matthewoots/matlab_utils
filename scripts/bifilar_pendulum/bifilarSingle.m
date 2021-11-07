clc
clear all
close all

g = 9.81;   % gravitational acceleration m/s^2
m = 0.1;    % mass of obj kg
d = 0.2;    % displacement between the 2 filars m
l = 0.2;    % length of the wire m
w = 0.2;    % angular frequency hz(1/s)

I = (m * g * d^2) / (4 * l * w^2);