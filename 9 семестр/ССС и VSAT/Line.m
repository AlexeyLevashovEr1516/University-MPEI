close all; clear all; clc;

x = 0:100;
y1 = x;
y2 = x/2;

figure(1);
hold on, grid on;
plot(x, y1, x, y2, 'LineWidth', 2)