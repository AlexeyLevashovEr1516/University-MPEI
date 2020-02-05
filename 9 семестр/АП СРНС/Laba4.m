close all; clear all; clc;

c = 3e8;
tay = 1e-6;

%R = [5.45; 9.99; 2.473];
R = [3.19; 2.423; 6.232];
X1 = [4.7; 5.09; 0.774];
X2 = [1.579; 0.858; 0.78];
X3 = [9; 3.476; 1.557];
%x = [0; 0; 0];
D = c*tay;
Ksi = [0; 0; 0; 0];

epsilon = 1e-3;
epsi = 10;
k_iter = 0;

while (epsi>=epsilon)
    f1 = Ksi(4) + sqrt((X1(1)-Ksi(1))^2 + (X1(2)-Ksi(2))^2 + (X1(3)-Ksi(3))^2);
    f2 = Ksi(4) + sqrt((X2(1)-Ksi(1))^2 + (X2(2)-Ksi(2))^2 + (X2(3)-Ksi(3))^2);
    f3 = Ksi(4) + sqrt((X3(1)-Ksi(1))^2 + (X3(2)-Ksi(2))^2 + (X3(3)-Ksi(3))^2);    
    f = [f1; f2; f3];
    
    %% первая строка
    H11 = (X1(1)-Ksi(1))/f1;
    H12 = (X1(2)-Ksi(2))/f1;
    H13 = (X1(3)-Ksi(3))/f1;
    %% вторая строка
    H21 = (X2(1)-Ksi(1))/f2;
    H22 = (X2(2)-Ksi(2))/f2;
    H23 = (X2(3)-Ksi(3))/f2;
    %% третья строка
    H31 = (X3(1)-Ksi(1))/f3;
    H32 = (X3(2)-Ksi(2))/f3;
    H33 = (X3(3)-Ksi(3))/f3;
    
    H = -[H11 H12 H13 1;
          H21 H22 H23 1;
          H31 H32 H33 1];
    
    Ksi_old = Ksi;
    Ksi = Ksi + inv(H'*H)*H'*(R - f);
    
    epsi = sqrt((Ksi(1)-Ksi_old(1))^2 + (Ksi(2)-Ksi_old(2))^2 + (Ksi(3)-Ksi_old(3))^2);   
    k_iter = k_iter + 1
end