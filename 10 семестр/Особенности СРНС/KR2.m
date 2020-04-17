close all; clear all; clc;
format long

Ft = 10.23e6; % ������� ������� �������� [���/�]
chip = 1/Ft % ������������ ������������� ������� [�]
L = 10230; % ����� ��������� ����� [���]
T = 1e-3; % ������ ��������� ����� [�]

j = 5; % ��������� ����� �� � ��

% ��� ������������ �� L3OCp ������������ ��1 � ��3
DA1 = [ 0 0 1 1 0 1 0 0 1 1 1 0 0 0 ]; % �������� ��������� ��1
DA3 = [ 1 0 0 0 1 0 1 ]; % �������� ��������� ��3

% ������������ ��
for k = 1:L
    % �����
    out01(k) = xor( DA1(14) , DA3(7) );
    
    % �������� �����
    fb_DA1 = xor( DA1(4),DA1(8) );
    fb_DA1 = xor( fb_DA1,DA1(13) );
    fb_DA1 = xor( fb_DA1,DA1(14) );
    fb_DA3 = xor( DA3(6),DA3(7) );
    
    % �����
    DA1(2:14) = DA1(1:13);
    DA3(2:7) = DA3(1:6);
    DA1(1) = fb_DA1;
    DA3(1) = fb_DA3;
end

first16b = out01(1:16)
last16b = out01(L-15:L)

% ��������� ������� +-1
for k = 1:L
    if out01(k)
        out(k) = -1;
    else
        out(k) = +1;
    end
end

% ������ ���
S = fft(out);
SS = S.*conj(S);
AKF = real( ifft(SS) );
AKF_plot = [AKF(L:-1:2),AKF];

% ������ ���������
Apeak = max(abs(AKF(2:L)))
Lpeak = 10*log10(Apeak/AKF(1))

Astd = std(AKF(2:L))
Lstd = 10*log10(Astd/AKF(1))

% �������
figure
hold on
grid on
plot([-L+1:L-1], AKF_plot)
xlim([-L L]);
ylim([-500 500]);
title('������������������ ������� R(\tau)')
ylabel('R(\tau)');
xlabel('\tau, �');