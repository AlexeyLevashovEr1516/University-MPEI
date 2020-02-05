function MeowMuuur()
    clear all;

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
    phi = -(1.25*pi);
    delta = -(1.5*pi);

    % Вектор информативных параметров

    lambda = [A1 A2 w phi delta];
    for i = [1:100]
        A1 = lambda(1);
        A2 = lambda(2);
        w = lambda(3);
        phi = lambda(4);
        delta = lambda(5);
        for k = [1:M-1]
            % Упрощение
            c1 = cos(w * k * Ts + phi);
            c2 = sin(w * k * Ts + phi);
            c3 = cos(w * k * Ts + phi + delta);
            c4 = sin(w * k * Ts + phi + delta);

            % Производные 1ого порядка функции правдоподбия:
            LnA1(k) = (1/sigma^2)*sum(c1*Y(:,1) + c2*Y(:,2)) - (1/sigma^2)*A1*M;
            LnA2 = (1/sigma^2)*sum(c3*Y(:,3) + c4*Y(:,4)) - (1/sigma^2)*A2*M;
            Lnw = (1/sigma^2)*sum(k*Ts*(-A1*c2*Y(:,1) + A1*c1*Y(:,2) - A2*c4*Y(:,3) + A2*c3*Y(:,4)));
            Ln_phi = (1/sigma^2)*sum(-A1*c2*Y(:,1) + A1*c1*Y(:,2) - A2*c4*Y(:,3) - A2*c3*Y(:,4));
            Ln_delta = (1/sigma^2)*sum(-A2*c4*Y(:,3) + A2*c3*Y(:,4));

            % Производные 2ого порядка функции правдоподбия:
            LnA1_2 = - (1/sigma^2)*M;
            LnA1A2 = 0;
            LnA2_2 = - (1/sigma^2)*M;
            LnA1w = (1/sigma^2)*sum(k*Ts*(-c2*Y(:,1) + c1*Y(:,2)));
            LnA1_phi = (1/sigma^2)*sum(-c2*Y(:,1) + c1*Y(:,2));
            LnA1_delta = 0;
            LnA2w = (1/sigma^2)*sum(k*Ts*(-c4*Y(:,3) + c3*Y(:,4)));
            LnA2_phi = (1/sigma^2)*sum(-c4*Y(:,3) + c3*Y(:,4));
            LnA2_delta = (1/sigma^2)*sum(-c4*Y(:,3) + c3*Y(:,4));
            Lnw_2 = (1/sigma^2)*sum(k^2*Ts^2*(-A1*c1*Y(:,1) - A1*c2*Y(:,2) - A2*c3*Y(:,3) - A2*c4*Y(:,4)));
            Lnw_phi = (1/sigma^2)*sum(k*Ts*(-A1*c1*Y(:,1) - A1*c2*Y(:,2) - A2*c3*Y(:,3) - A2*c4*Y(:,4)));
            Lnw_delta = (1/sigma^2)*sum(k*Ts*(-A2*c4*Y(:,3) + A2*c3*Y(:,4)));
            Ln_phi_2 = (1/sigma^2)*sum(-A1*c1*Y(:,1) - A1*c2*Y(:,2) - A2*c3*Y(:,3) - A2*c4*Y(:,4));
            Ln_delta_2 = (1/sigma^2)*sum(-A2*c3*Y(:,3) - A2*c4*Y(:,4));
            Ln_phi_delta = (1/sigma^2)*sum(-A2*c3*Y(:,3) - A2*c4*Y(:,4));
        end
        
        Ln = [sum(LnA1) LnA2 Lnw Ln_phi Ln_delta];
        Ln_2 = [LnA1_2 LnA1A2 LnA1w LnA1_phi LnA1_delta; LnA1A2 LnA2_2 LnA2w LnA2_phi LnA2_delta; LnA1w LnA2w Lnw_2 Lnw_phi Lnw_delta; ...
                LnA1_phi LnA2_phi Lnw_phi Ln_phi_2 Ln_phi_delta; LnA1_delta LnA2_delta Lnw_delta Ln_phi_delta Ln_delta_2];
        lambda = lambda - Ln * inv(Ln_2)
        i
    end
end