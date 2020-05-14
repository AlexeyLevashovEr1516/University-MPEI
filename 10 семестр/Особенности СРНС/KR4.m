close all; clear all; clc;
format long

delta_f = 2.046e6; % [√ц] полоса приемника


f = -delta_f/2:1:delta_f/2; % [√ц] полоса частот

tau_L2CM = 1e-6/0.5115;
S_L2CM = tau_L2CM / 2 * sinc(pi*f*tau_L2CM).^2;

tau_L2CL = 1e-6/0.5115;
S_L2CL = tau_L2CL / 2 * sinc(pi*f*tau_L2CL).^2;

tau_L2PY = 1e-6/10.23;
S_L2PY = tau_L2PY / 2 * sinc(pi*f*tau_L2PY).^2;

tau_L2M = 1e-6/5.115;
fs_L2M = 10*1.023e6;
S_L2M = tau_L2M / 2 * sinc(pi*f*tau_L2M).^2 .* tan(pi/2 * f/fs_L2M).^2;

f = f*1e-6;
figure
hold on
grid on
plot(f, S_L2CM, f, S_L2CL, f, S_L2PY, f, S_L2M)
xlim([f(1) f(end)]);
xlabel('f, MHz');
ylabel('S(f)');
legend('S(f) L2CM','S(f) L2CL','S(f) L2P(Y)','S(f) L2M');

Ksd1 = sum(S_L2M.*S_L2CL);
Ksd2 = sum(S_L2PY.*S_L2CL);
Ksd3 = sum(S_L2CM.*S_L2CL);
Ksd4 = sum(S_L2CL.*S_L2CL);

q_JN0 = 45;
q_JN0 = 10^(0.1*q_JN0);

Kjam = 1 / (1 + q_JN0*((10*Ksd1) + (16*Ksd2) + (10*Ksd3) + (9*Ksd4))) 



