close all; clear all; clc;
format long

% GPS L2C
f_L2 = 1227.6e6; % Несущая частота [Гц]

% CM-код (информационная компонента)
T_cm = 20e-3; % Период кода [с]
L_cm = 10230; % Длина кода [бит]
Ft_cm = 0.5115e6; % Частота выборки символов [бит/с]
chip_cm = 1/Ft_cm; % Длительность элементарного символа [с]

% CL-код (пилотная компонента)
T_cl = 1.5; % Период кода [с]
L_cl = 767250; % Длина кода [бит]
Ft_cl = 0.5115e6; % Частота выборки символов [бит/с]
chip_cl = 1/Ft_cl; % Длительность элементарного символа [с]

% Нав. сообщение
chip_nd = 20e-3; % Длительность элементарного символа [с]
G_nav_data = [-1 1]; % Содержание нав. сообщения

% Формирования ДК CM-код
reg_cm = -ones(1,27); % Начально состояние регистра
for k = 1:L_cm
    % Выход
    G_cm(k) = reg_cm(27);
    
    % Обратная связь
    feedback = reg_cm(3)*reg_cm(6)*reg_cm(8)*reg_cm(11)*reg_cm(14)*reg_cm(16)*reg_cm(18)*reg_cm(21)*reg_cm(22)*reg_cm(23)*reg_cm(24)*reg_cm(27);
    
    % Сдвиг
    reg_cm(2:27) = reg_cm(1:26);
    reg_cm(1) = feedback;
end

% Формирования ДК CL-код
reg_cl = -ones(1,27); % Начально состояние регистра
for k = 1:L_cl
    % Выход
    G_cl(k) = reg_cl(27);
    
    % Обратная связь
    feedback = reg_cl(3)*reg_cl(6)*reg_cl(8)*reg_cl(11)*reg_cl(14)*reg_cl(16)*reg_cl(18)*reg_cl(21)*reg_cl(22)*reg_cl(23)*reg_cl(24)*reg_cl(27);
    
    % Сдвиг
    reg_cl(2:27) = reg_cl(1:26);
    reg_cl(1) = feedback;
end

% Натягивание на время
delta_f = 2.046e6;
fs = 4 * delta_f;
ts = 1 / fs;
f0 = fs / 4;
A = 1;

tos = 0;
tof = 20e-3;
t = tos:ts:tof;
N = length(t);

k0_cm = 1;
k1_cm = 0;
G1_cm = G_cm(1);

for k = 1:N
    tk = ts*k;
    chipk_cm = k0_cm * chip_cm;
    
    if tk == chipk_cm
        k0_cm = k0_cm + 1;
        if k0_cm == L_cm*(k1_cm+1)
            k1_cm = k1_cm + 1;
        end
        if k1_cm == 0
            Gk = k0_cm;
        else
            Gk = k0_cm - k1_cm*L_cm + 1;
        end
        G1_cm = G_cm(Gk);
    end
    Gcm(k) = G1_cm;
end

k0_cl = 1;
k1_cl = 0;
G1_cl = G_cl(1);

for k = 1:N
    tk = ts*k;
    chipk_cl = k0_cl * chip_cm;
    
    if tk == chipk_cl
        k0_cl = k0_cl + 1;
        if k0_cl == L_cl*(k1_cl+1)
            k1_cl = k1_cl + 1;
        end
        if k1_cl == 0
            Gk = k0_cl;
        else
            Gk = k0_cl - k1_cl*L_cl + 1;
        end
        G1_cl = G_cl(Gk);
    end
    Gcl(k) = G1_cl;
end

s_L2C = A*Gcm.*cos(2*pi*f0*t);

% Графики
% figure
% subplot(3,1,1)
% hold on
% grid on
% plot(t, Gcm)
% title('ДК CM')
% 
% subplot(3,1,2)
% hold on
% grid on
% plot(t, Gcl)
% title('ДК CL')
% 
% subplot(3,1,3)
% hold on
% grid on
% plot(t, s_L2C)
% title('Сигнал L2C')

figure
hold on
grid on
plot(t, s_L2C, ':')
plot(t, Gcm)