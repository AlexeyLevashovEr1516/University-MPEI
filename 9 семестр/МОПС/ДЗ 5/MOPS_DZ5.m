close all; clear all; clc;
format long

%% Параметры
% От НАП:
T = 10e-3;
Td = 0.2e-6;
N = T/Td;

t_start = 0;
t_stop = 10;
t = t_start:T:t_stop;
td = t_start:Td:t_stop;

q_dB = 30;
q = 10^(q_dB/10);

alpha = 1;
sigma_alpha = 10;
c = 3e8;

f0 = 1602e6;
omega0 = 2*pi*f0;
fp = 2e3;
omegap = 2*pi*fp;

% От ИНС
alpha_delta = 0.1;
sigma_delta = 1;

gamma = 0;

%% Шум от ИНС:
S_chi = 2*sigma_delta^2*alpha_delta;
D_chi = S_chi/(2*T);

%% Шум наблюдений
a0 = 1;
sigma_n = a0/(2*sqrt(q*Td));
D_n = sigma_n^2;

%% Формирующий шум
S_xi = 2*sigma_alpha^2*alpha*(omega0/c)^2;
D_xi = S_xi/(2*T);
D_zeta = 0.5^2;

%% Коэффициенты фильтра
F = [1 0 0 0;
     0 1 T 0;
     0 0 1 -(omega0/c)*T;
     0 0 0 1-alpha_delta*T];
G = [T 0;
     0 0;
     0 0;
     0 alpha_delta*T];
C = [1 0 0 0;
     0 1 0 0];
D_f = [D_zeta 0;
       0 D_chi];

%% Начальные условия
x = [1; pi/12; 100; 1];
D = [.3^2 0 0 0;
     0 pi^2 0 0;
     0 0 34^2 0;
     0 0 0 1^2];
xf = [.5; 0; 0; 0];

%% Выделение памяти и начальные приближения
a = nan(size(t));           a(1) = x(1);
phi = nan(size(t));         phi(1) = x(2);
OMEGA = nan(size(t));       OMEGA(1) = x(3);
a_extr = nan(size(t));      a_extr(1) = 0;
phi_extr = nan(size(t));    phi_extr(1) = 0;
OMEGA_extr = nan(size(t));  OMEGA_extr(1) = 0;
D22 = nan(size(t));         D22(1) = D(2,2);

for k = 2:length(t)   
    x = F*x + G*randn(1,1)*sqrt([D_zeta; D_chi]); %% Марковский СП
    x(3) = x(3) + gamma*T;
    a(k)     = x(1);
    phi(k)   = x(2);
    OMEGA(k) = x(3);
    delta(k) = x(4);
    
    v_true(k) = gamma*(c/omega0) - x(4);
    gamma = gamma + randn(1,1)*sqrt(0.0001*k*T); %% Винеровский СП
    gammak(k) = gamma;
    
    %% экстраполяция
    xf = F*xf;
    D = F*D*F' + G*D_f*G';		
    W = N/(2*D_n)*[1 0;
                   0 xf(1)^2];
               
    a_extr(k) = xf(1);
    phi_extr(k) = xf(2);
    OMEGA_extr(k) = xf(3);
               
    %% Дискриминация
    for i = 1:N
        i_m = (k-2)*N+i;
        if i_m*Td <= 5
            ai = 1;
        else
            ai = 0.5;
        end
        y = ai*cos(omegap*i_m*Td + phi(k)) + randn(1,1)*sigma_n;
        I(i) = y*cos(omegap*i_m*Td + phi_extr(k));
        Q(i) = y*sin(omegap*i_m*Td + phi_extr(k));
    end
    U_d1 = sum(I) * (1/D_n) - (xf(1)*N)/(2*D_n);
    U_d2 = -sum(Q) * (xf(1)/D_n);
    u_d = [U_d1; U_d2];

    %% Оценка 
    D = inv(inv(D) + C'*W*C);
    xf = xf + D*C'*u_d;

    D22(k) = D(2,2);
end

epsilon_phi = (phi_extr - phi);

figure(1);
hold on, grid on;
plot(t, epsilon_phi*180/pi, t, [3*sqrt(D22); -3*sqrt(D22)]*180/pi)
ylim([-100 100]);
xlabel("Время, с");
ylabel("Мгновенные ошибки фильтрации фазы, град");

figure(2);
hold on, grid on;
plot(t, v_true, t, delta)
xlabel("Время, с");
ylabel("Реализация истинного ускорения и погрешность измерений радиального ускорения");

figure(3);
hold on, grid on;
plot(t, gammak)
xlabel("Время, с");
ylabel("Реализация истинного ускорения и погрешность измерений радиального ускорения");