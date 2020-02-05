clear all; clc; close all;

%I = [sqrt(1/2); sqrt(1/2); -sqrt(1/2); -sqrt(1/2); 1; -1; 0; 0];
%Q = [sqrt(1/2); -sqrt(1/2); sqrt(1/2); -sqrt(1/2); 0; 0; 1; -1];

I = [sqrt(1/2); sqrt(1/2); -sqrt(1/2); -sqrt(1/2); 1; -1; 0; 0];
Q = [sqrt(1/2); -sqrt(1/2); sqrt(1/2); -sqrt(1/2); 0; 0; 1; -1];

maxIQ = 1.4*max(abs(I + 1i*Q));

figure(1);
plot(I, Q, '*')
hold on
plot([-maxIQ maxIQ], [0 0], 'k');
plot([0 0], [-maxIQ maxIQ], 'k');
quiver(-maxIQ, 0, 2*maxIQ, 0, 1, 'k');
quiver(0, -maxIQ, 0, 2*maxIQ, 1, 'k');
hold off
xlabel('I')
ylabel('Q');
axis equal