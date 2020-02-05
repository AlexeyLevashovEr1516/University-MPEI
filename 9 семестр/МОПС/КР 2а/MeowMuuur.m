function MeowMuuur()
    clear all;
    format long
    % Построение графика из выборки в текстовом файле
    Y=load('Input_Y0toT.txt');

    % Исходные данные
    M = 2048;
    sigma = 10;
    Fs = 47.5e6;
    Ts = 1/Fs;
    k = 1:M-1;
    t = Ts * (0:M-1);

    % Задание начальных приближений
    A1 = 6100;
    A2 = 6800;
    f = 2.5e5;
    w = 2*pi*f;
%     phi = -(1.25*pi);
    phi = 2.319659512960347;
   % delta = -(1.5*pi);
   delta = 1.546440257083316;

    % Вектор информативных параметров

    lambda = [A1 A2 w phi delta];
    for i = [1:100]
        A1 = lambda(1);
        A2 = lambda(2);
        w = lambda(3);
        phi = lambda(4);
        delta = lambda(5);
        for k = [1:M]
            % Упрощение
            c1 = cos(w * k * Ts + phi);
            c2 = sin(w * k * Ts + phi);
            c3 = cos(w * k * Ts + phi + delta);
            c4 = sin(w * k * Ts + phi + delta);

            % Производные 1ого порядка функции правдоподбия:
            LnA1(k) = (1/sigma^2)*(c1*Y(k,1) + c2*Y(k,2)) - (1/sigma^2)*A1;
            LnA2(k) = (1/sigma^2)*(c3*Y(k,3) + c4*Y(k,4)) - (1/sigma^2)*A2;
            Lnw(k) = (1/sigma^2)*(k*Ts*(-A1*c2*Y(k,1) + A1*c1*Y(k,2) - A2*c4*Y(k,3) + A2*c3*Y(k,4)));
            Ln_phi(k) = (1/sigma^2)*(-A1*c2*Y(k,1) + A1*c1*Y(k,2) - A2*c4*Y(k,3) - A2*c3*Y(k,4));
            Ln_delta(k) = (1/sigma^2)*(-A2*c4*Y(k,3) + A2*c3*Y(k,4));

            % Производные 2ого порядка функции правдоподбия:
            LnA1_2(k) = - (1/sigma^2);
            LnA1A2(k) = 0;
            LnA2_2(k) = - (1/sigma^2);
            LnA1w(k) = (1/sigma^2)*(k*Ts*(-c2*Y(k,1) + c1*Y(k,2)));
            LnA1_phi(k) = (1/sigma^2)*(-c2*Y(k,1) + c1*Y(k,2));
            LnA1_delta(k) = 0;
            LnA2w(k) = (1/sigma^2)*(k*Ts*(-c4*Y(k,3) + c3*Y(k,4)));
            LnA2_phi(k) = (1/sigma^2)*(-c4*Y(k,3) + c3*Y(k,4));
            LnA2_delta(k) = (1/sigma^2)*(-c4*Y(k,3) + c3*Y(k,4));
            Lnw_2(k) = (1/sigma^2)*(k^2*Ts^2*(-A1*c1*Y(k,1) - A1*c2*Y(k,2) - A2*c3*Y(k,3) - A2*c4*Y(k,4)));
            Lnw_phi(k) = (1/sigma^2)*(k*Ts*(-A1*c1*Y(k,1) - A1*c2*Y(k,2) - A2*c3*Y(k,3) - A2*c4*Y(k,4)));
            Lnw_delta(k) = (1/sigma^2)*(k*Ts*(-A2*c4*Y(k,3) + A2*c3*Y(k,4)));
            Ln_phi_2(k) = (1/sigma^2)*(-A1*c1*Y(k,1) - A1*c2*Y(k,2) - A2*c3*Y(k,3) - A2*c4*Y(k,4));
            Ln_delta_2(k) = (1/sigma^2)*(-A2*c3*Y(k,3) - A2*c4*Y(k,4));
            Ln_phi_delta(k) = (1/sigma^2)*(-A2*c3*Y(k,3) - A2*c4*Y(k,4));
        end
        
        Ln = [sum(LnA1) sum(LnA2) sum(Lnw) sum(Ln_phi) sum(Ln_delta)];
        Ln_2 = [sum(LnA1_2) sum(LnA1A2) sum(LnA1w) sum(LnA1_phi) sum(LnA1_delta);
                sum(LnA1A2) sum(LnA2_2) sum(LnA2w) sum(LnA2_phi) sum(LnA2_delta);
                sum(LnA1w) sum(LnA2w) sum(Lnw_2) sum(Lnw_phi) sum(Lnw_delta); ...
                sum(LnA1_phi) sum(LnA2_phi) sum(Lnw_phi) sum(Ln_phi_2) sum(Ln_phi_delta);
                sum(LnA1_delta) sum(LnA2_delta) sum(Lnw_delta) sum(Ln_phi_delta) sum(Ln_delta_2)];
        lambda = lambda - Ln * inv(Ln_2)
        i
    end
end