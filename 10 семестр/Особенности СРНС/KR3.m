close all; clear all; clc;
format long

%% GPS L2C
f_L2 = 1227.6e6; % Несущая частота [Гц]

%% Нав. сообщение Navigation Data
chip_ND = 20e-3; % Длительность элементарного символа [с]
G_ND = [-1 1]; % Содержание нав. сообщения
L_ND = length(G_ND);

%% Формирование ДК 
% CM-код (информационная компонента)
T_CM = 20e-3; % Период кода [с]
L_CM = 10230; % Длина кода [бит]
Ft_CM = 0.5115e6; % Частота выборки символов [бит/с]
chip_CM = 1/Ft_CM; % Длительность элементарного символа [с]

% CL-код (пилотная компонента)
T_CL = 1.5; % Период кода [с]
L_CL = 767250; % Длина кода [бит]
Ft_CL = 0.5115e6; % Частота выборки символов [бит/с]
chip_CL = 1/Ft_CL; % Длительность элементарного символа [с]

% Начальное состояние регистров сдвига (Initial Shift Register State)
ISRS_CM = [1 1 1 1 0 1 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 1]; %756014035
ISRS_CL = [1 0 1 0 0 0 1 1 0 1 1 0 0 0 1 0 0 0 0 1 1 1 1 0 0 1 0]; %506610362

% Расчет ДК
DK_L2C_CM_out = DK_L2C_calc( ISRS_CM,L_CM );
DK_L2C_CL_out = DK_L2C_calc( ISRS_CL,L_CL );

%% Натягивание на время
fc = 2.046e6; % Ширина спектра по главномым лепесткам - частота выборки символов ДК
fs = 4 * fc; % Частота дискретизации
ts = 1 / fs; % Период дискретизации
f0 = fs / 4; % Промежуточная частота 
A = 1; % Амплитуда

k = 0; % номер текущего отсчета
tk = 0; % Время старта
toe = 20e-3; % Время окончания 
t_akf = [-toe+ts:ts:toe-ts];

while tk <= toe
    k = k + 1;
    tout(k) = tk;
    
    %Формирование CM-кода
    N_chip_CM(k) = mod( fix(tk/chip_CM), L_CM ) + 1;
    DK_out_CM(k) = DK_L2C_CM_out(N_chip_CM(k));
    
    % Формирование CL-кода
    N_chip_CL(k) = mod( fix(tk/chip_CL), L_CL ) + 1;
    DK_out_CL(k) = DK_L2C_CL_out(N_chip_CL(k));
    
    % Формирование навигационного сообщения
    N_chip_ND(k) = mod( fix(tk/chip_ND), L_ND ) + 1;
    DK_out_ND(k) = G_ND(N_chip_ND(k));
    
    % временное уплотнение Time Multiplexing
    TM_valid(k) = mod( fix(2*(tk/chip_CM)), 2 );
    if TM_valid(k)
        DKout(k) = DK_out_CL(k);
    else
        DKout(k) = DK_out_CM(k)*DK_out_ND(k);
    end

    tk = tk + ts;
end

signal_L2C = A*DKout.*cos(2*pi*f0*tout);

%% Энергетический спектр и АКФ
% ДК CL
S_CL = fft(DK_out_CL);
SS_CL = S_CL.*conj(S_CL);
AKF_CL = real( ifft(SS_CL) );
AKF_CL_plot = [AKF_CL(length(AKF_CL):-1:2),AKF_CL];

% ДК CM
S_CM = fft(DK_out_CM);
SS_CM = S_CM.*conj(S_CM);
AKF_CM = real( ifft(SS_CM) );
AKF_CM_plot = [AKF_CM(length(AKF_CM):-1:2),AKF_CM];

% Сигнал
F = 0:1/toe:(fs-1/toe); % Формирование оси частот

S_signal = fft(signal_L2C);
SS_signal = S_signal.*conj(S_signal);
SSS_signal = 2*SS_signal(1:length(F));

AKF_signal = real( ifft(SSS_signal) );
AKF_signal_plot = [AKF_signal(length(AKF_signal):-1:2),AKF_signal];

%% Графики
figure
hold on
grid on
plot(tout*1e6, signal_L2C)
xlim([0 5]*chip_CM*1e6)
title('Сигнал L2C')
xlabel('Время, мкс')
ylabel('Амплитуда, В')

figure
hold on
grid on
plot(tout*1e3,DKout)
xlabel('Время, мс')

figure
hold on
grid on
plot(tout*1e3,N_chip_CM,tout*1e3,N_chip_CL,tout*1e3,TM_valid)
legend('DK CM','DK CL','TM valid')
xlabel('Время, мс')

figure
hold on
grid on
plot(F*1e-6, 10*log10(abs(SSS_signal)/abs(max(SSS_signal))))
xlim([0 2*fc]*1e-6)
title('Энергеический спектр сигнала GPS L2C')
xlabel('Частота, МГц')
ylabel('S(f), дБ')

figure
hold on
grid on
plot(t_akf*1e3, AKF_signal_plot)
title('Автокорреляционная функция R(\tau) сигнала GPS L2C')
ylabel('R(\tau)');
xlabel('\tau, мс');
% 
% figure
% hold on
% grid on
% plot([-length(AKF_CM)+1:length(AKF_CM)-1], AKF_CM_plot)
% xlim([-length(AKF_CM) length(AKF_CM)]);
% title('Автокорреляционная функция R(\tau) CM кода')
% ylabel('R(\tau)');
% xlabel('\tau, с');
% 
% figure
% hold on
% grid on
% plot([-length(AKF_CL)+1:length(AKF_CL)-1], AKF_CL_plot)
% xlim([-length(AKF_CL) length(AKF_CL)]);
% title('Автокорреляционная функция R(\tau) CL кода')
% ylabel('R(\tau)');
% xlabel('\tau, с');