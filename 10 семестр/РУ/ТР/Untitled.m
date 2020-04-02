clear all; close all; tic;
format long;

A = 15;
B = 7;
Tf = 0.1;
a = 0.4;
Ky2 = 3.5;
T1 = 0.5;
alpha = 0.5;
Tym = 0.01;
Kgp = 1;
Sd = 6;
%Sksi = 1e-4;

K0(1) = Sd;
K1(1) = Sd;

Sksi = 0:0.001:20;
N = length(Sksi);

for i = 1:N
   m_theta(i) = alpha / (K0(i) * Ky2 * Kgp);
   %D_theta(i) = Sksi(i) * ((-T1^2 * Kgp + (1 / (K1(i) * Ky2)))/(2 * K1(i) * Ky2 * (Tym - T1)));
   D_theta(i) = Sksi(i) * ((T1^2 * K1(i) * Ky2 * Kgp + 1) / (2 * K1(i)^2 * (T1 - Tym)));
   sigma_theta(i) = sqrt(D_theta(i));
   
   K0(i+1) = (A / m_theta(i)) * sin(alpha * m_theta(i)) * exp(-((alpha * sigma_theta(i))^2)/2);
   K1(i+1) = A * alpha * cos(alpha * m_theta(i)) * exp(-((alpha * sigma_theta(i))^2)/2);
end

theta = m_theta + 3*sigma_theta;

figure
subplot(3,1,1)
hold on; grid on;
ylim([0 0.04]);
plot(Sksi,m_theta);
title('Зависимость m_\theta(S_\xi)')
ylabel('m_\theta, град');
xlabel('S_\xi, В^2с');

subplot(3,1,2)
hold on; grid on;
ylim([0 2.5]);
plot(Sksi,sigma_theta);
title('Зависимость \sigma_\theta(S_\xi)')
ylabel('\sigma_\theta, град');
xlabel('S_\xi, В^2с');

subplot(3,1,3)
hold on; grid on;
ylim([0 7]);
plot(Sksi,theta);
title('Зависимость \theta_{max}(S_\xi)')
ylabel('\theta_{max}, град');
xlabel('S_\xi, В^2с');
