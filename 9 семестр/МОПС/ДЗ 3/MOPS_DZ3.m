function Main()
    close all; clear all; clc;
    
    %% Параметры
    T = 10e-3;
    t = 0:T:10;
    q_dB = 30;
    q = 10^(q_dB/10);
    alpha = 1;
    c = 3e8;
    f0 = 1602e6;
    omega0 = 2*pi*f0;
    sigma_alpha = 10;
    
    %% СПМ шума наблюдений
    N0 = 2/(q*T^2)*(1+1/(2*q*T));
    D_n = N0/(2*T);
    
    %% СПМ формирующего шума
    S = 2*sigma_alpha^2*alpha*(omega0/c)^2;
    D_ksi = S/(2*T);
    
    %% Коэффициенты фильтра
    H = [1 0];
    F = [1 T; 0 (1-alpha*T)];
    G = [0; alpha*T];
    
    %% Начальное приближение
    x = [100; 100];
    D0 = [34^2 0; 0 340^2];
    xf = [0; 0];
          
    %% Выделение памяти и начальные приближения
    OMEGA       = nan(size(t)); OMEGA(1)      = x(1);
    Nu          = nan(size(t)); Nu(1)         = x(2);
    D_OMEGA     = nan(size(t)); D_OMEGA(1)    = D0(1,1);
    y           = nan(size(t)); y(1)          = 0;
    OMEGA_extr  = nan(size(t)); OMEGA_extr(1) = xf(1);
    
    for k = 2:length(t)
        ksi = randn(1,1)*sqrt(D_ksi);
        x = F*x + G*ksi;
        OMEGA(k) = x(1);
        Nu(k)    = x(2);
        
        %% экстраполяция
        D0 = F*D0*F' + G*D_ksi*G';
        K = D0*H'*inv(H*D0*H' + D_n);
        xf = F * xf;
        
        %% наблюдения
        yk = OMEGA(k) + randn(1,1)*sqrt(D_n);
        y(k) = yk;
        
        %% коррекция
        D0 = (eye(length(x)) - K*H)*D0;
        D_OMEGA(k) = D0(1,1);
        xf = xf + K * (yk - H*xf);
        OMEGA_extr(k) = xf(1);
        
    end
    
    d_OMEGA = (OMEGA_extr - OMEGA);
    std(d_OMEGA)^2 % Для всего временного участка
    mean(D_OMEGA())
    std(d_OMEGA(30:1000))^2 % Для установившегося режима
    mean(D_OMEGA(30:1000))
    
    figure(1);
    hold on, grid on;
    plot(t, OMEGA)
    xlabel("Время, с");
    ylabel("Доплеровская частота, Гц");
    
    figure(2);
    hold on, grid on;
    plot(t, sqrt(D_OMEGA))
    xlabel("Время, с");
    ylabel("СКО фильтрации частоты");
    
    figure(3);
    hold on, grid on;
    plot(t, d_OMEGA, t, [+3*sqrt(D_OMEGA); -3*sqrt(D_OMEGA)])
    xlabel("Время, с");
    ylabel("Мгновенная ошибка фильтрации частоты, Гц");

end