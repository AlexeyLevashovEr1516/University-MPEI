close all; clear all; clc;
format long

%% GPS L2C
f_L2 = 1227.6e6; % ������� ������� [��]

%% ���. ��������� Navigation Data
chip_ND = 20e-3; % ������������ ������������� ������� [�]
G_ND = [-1 1]; % ���������� ���. ���������
L_ND = length(G_ND);

%% ������������ �� 
% CM-��� (�������������� ����������)
T_CM = 20e-3; % ������ ���� [�]
L_CM = 10230; % ����� ���� [���]
Ft_CM = 0.5115e6; % ������� ������� �������� [���/�]
chip_CM = 1/Ft_CM; % ������������ ������������� ������� [�]

% CL-��� (�������� ����������)
T_CL = 1.5; % ������ ���� [�]
L_CL = 767250; % ����� ���� [���]
Ft_CL = 0.5115e6; % ������� ������� �������� [���/�]
chip_CL = 1/Ft_CL; % ������������ ������������� ������� [�]

% ��������� ��������� ��������� ������ (Initial Shift Register State)
ISRS_CM = [1 1 1 1 0 1 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 1]; %756014035
ISRS_CL = [1 0 1 0 0 0 1 1 0 1 1 0 0 0 1 0 0 0 0 1 1 1 1 0 0 1 0]; %506610362

% ������ ��
DK_L2C_CM_out = DK_L2C_calc( ISRS_CM,L_CM );
DK_L2C_CL_out = DK_L2C_calc( ISRS_CL,L_CL );

%% ����������� �� �����
fc = 2.046e6; % ������ ������� �� ��������� ��������� - ������� ������� �������� ��
fs = 4 * fc; % ������� �������������
ts = 1 / fs; % ������ �������������
f0 = fs / 4; % ������������� ������� 
A = 1; % ���������

k = 0; % ����� �������� �������
tk = 0; % ����� ������
toe = 20e-3; % ����� ��������� 
t_akf = [-toe+ts:ts:toe-ts];

while tk <= toe
    k = k + 1;
    tout(k) = tk;
    
    %������������ CM-����
    N_chip_CM(k) = mod( fix(tk/chip_CM), L_CM ) + 1;
    DK_out_CM(k) = DK_L2C_CM_out(N_chip_CM(k));
    
    % ������������ CL-����
    N_chip_CL(k) = mod( fix(tk/chip_CL), L_CL ) + 1;
    DK_out_CL(k) = DK_L2C_CL_out(N_chip_CL(k));
    
    % ������������ �������������� ���������
    N_chip_ND(k) = mod( fix(tk/chip_ND), L_ND ) + 1;
    DK_out_ND(k) = G_ND(N_chip_ND(k));
    
    % ��������� ���������� Time Multiplexing
    TM_valid(k) = mod( fix(2*(tk/chip_CM)), 2 );
    if TM_valid(k)
        DKout(k) = DK_out_CL(k);
    else
        DKout(k) = DK_out_CM(k)*DK_out_ND(k);
    end

    tk = tk + ts;
end

signal_L2C = A*DKout.*cos(2*pi*f0*tout);

%% �������������� ������ � ���
% �� CL
S_CL = fft(DK_out_CL);
SS_CL = S_CL.*conj(S_CL);
AKF_CL = real( ifft(SS_CL) );
AKF_CL_plot = [AKF_CL(length(AKF_CL):-1:2),AKF_CL];

% �� CM
S_CM = fft(DK_out_CM);
SS_CM = S_CM.*conj(S_CM);
AKF_CM = real( ifft(SS_CM) );
AKF_CM_plot = [AKF_CM(length(AKF_CM):-1:2),AKF_CM];

% ������
F = 0:1/toe:(fs-1/toe); % ������������ ��� ������

S_signal = fft(signal_L2C);
SS_signal = S_signal.*conj(S_signal);
SSS_signal = 2*SS_signal(1:length(F));

AKF_signal = real( ifft(SSS_signal) );
AKF_signal_plot = [AKF_signal(length(AKF_signal):-1:2),AKF_signal];

%% �������
figure
hold on
grid on
plot(tout*1e6, signal_L2C)
xlim([0 5]*chip_CM*1e6)
title('������ L2C')
xlabel('�����, ���')
ylabel('���������, �')

figure
hold on
grid on
plot(tout*1e3,DKout)
xlabel('�����, ��')

figure
hold on
grid on
plot(tout*1e3,N_chip_CM,tout*1e3,N_chip_CL,tout*1e3,TM_valid)
legend('DK CM','DK CL','TM valid')
xlabel('�����, ��')

figure
hold on
grid on
plot(F*1e-6, 10*log10(abs(SSS_signal)/abs(max(SSS_signal))))
xlim([0 2*fc]*1e-6)
title('������������� ������ ������� GPS L2C')
xlabel('�������, ���')
ylabel('S(f), ��')

figure
hold on
grid on
plot(t_akf*1e3, AKF_signal_plot)
title('������������������ ������� R(\tau) ������� GPS L2C')
ylabel('R(\tau)');
xlabel('\tau, ��');
% 
% figure
% hold on
% grid on
% plot([-length(AKF_CM)+1:length(AKF_CM)-1], AKF_CM_plot)
% xlim([-length(AKF_CM) length(AKF_CM)]);
% title('������������������ ������� R(\tau) CM ����')
% ylabel('R(\tau)');
% xlabel('\tau, �');
% 
% figure
% hold on
% grid on
% plot([-length(AKF_CL)+1:length(AKF_CL)-1], AKF_CL_plot)
% xlim([-length(AKF_CL) length(AKF_CL)]);
% title('������������������ ������� R(\tau) CL ����')
% ylabel('R(\tau)');
% xlabel('\tau, �');