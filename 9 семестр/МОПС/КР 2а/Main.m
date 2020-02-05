function Main()
    close all; clear all;
    
    %% Задание параметров
    M = 2048;
    fs = 47.5e6;
    ts = 1 / fs;
    k = (0:M-1).';
    t = ts*(0:M-1).';
    n = 4;
    sigma = 10;
    
    %% Формирование наблюдения
    file = fopen('Input_Y0toT.txt');
    Y = [];
    while (~feof(file))
        scan = fscanf (file, '%f %f %f %f', [4 M]);
        Y = [Y; scan];
    end
    fclose(file);
    Y = Y';
    for i = 1:M
        Y1(1,i) = Y(i,1);
        Y2(1,i) = Y(i,2);
        Y3(1,i) = Y(i,3);
        Y4(1,i) = Y(i,4);
    end
    
    %% Начальное приближение
    A1 = 6100;
    A2 = 6800;
    f = 0.25e6;
    omega = 2 * pi * f;
    phi0 = - (1.25 * pi);
    delta_phi = - (1.5 * pi);
    lambda = [A1 A2 omega phi0 delta_phi]; 
    lambda_fist = lambda;
    
    %% Отношение правдоподобия
    for i = 2:12
    %while ((lambda_old - lambda) > (1e-5))
        for k = 1:M
        % Первые производные по параметрам
        d_lambda_1(k) = (1 / (sigma ^ 2)) * (cos(lambda(3) * k * ts + lambda(4)) * Y1(k) + ...
                                             sin(lambda(3) * k * ts + lambda(4)) * Y2(k));
        d_lambda_2(k) = (1 / (sigma ^ 2)) * (cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k) + ...
                                             sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));
        d_lambda_3(k) = (1 / (sigma ^ 2)) * (- lambda(1) * k * ts * sin(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                             + lambda(1) * k * ts * cos(lambda(3) * k * ts + lambda(4)) * Y2(k)...
                                             - lambda(2) * k * ts * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                             + lambda(2) * k * ts * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));
        d_lambda_4(k) = (1 / (sigma ^ 2)) * (- lambda(1) * sin(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                             + lambda(1) * cos(lambda(3) * k * ts + lambda(4)) * Y2(k)...
                                             - lambda(2) * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                             + lambda(2) * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));       
        d_lambda_5(k) = (1 / (sigma ^ 2)) * (- lambda(2) * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                             + lambda(2) * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));   
        % Вторые производные по параметрам
        % Первый столбец матрицы
        dd_lambda_11(k) = (1 / (sigma ^ 2)) * (-1);
        dd_lambda_21(k) = 0; 
        dd_lambda_31(k) = (1 / (sigma ^ 2)) * (- k * ts * sin(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                               + k * ts * cos(lambda(3) * k * ts + lambda(4)) * Y2(k));
        dd_lambda_41(k) = (1 / (sigma ^ 2)) * (- sin(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                               + cos(lambda(3) * k * ts + lambda(4)) * Y2(k));
        dd_lambda_51(k) = 0;
        
        % Второй столбец матрицы
        dd_lambda_12(k) = dd_lambda_21(k);                        
        dd_lambda_22(k) = (1 / (sigma ^ 2)) * (-1);
        dd_lambda_32(k) = (1 / (sigma ^ 2)) * (- k * ts * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               + k * ts * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));        
        dd_lambda_42(k) = (1 / (sigma ^ 2)) * (- sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               + cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k)); 
        dd_lambda_52(k) = (1 / (sigma ^ 2)) * (- sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               + cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k)); 

        % Третий столбец матрицы
        dd_lambda_13(k) = dd_lambda_31(k);                        
        dd_lambda_23(k) = dd_lambda_32(k);                      
        dd_lambda_33(k) = (1 / (sigma ^ 2)) * (- lambda(1) * ((k * ts) ^ 2) * cos(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                               - lambda(1) * ((k * ts) ^ 2) * sin(lambda(3) * k * ts + lambda(4)) * Y2(k)...
                                               - lambda(2) * ((k * ts) ^ 2) * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               - lambda(2) * ((k * ts) ^ 2) * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));        
        dd_lambda_43(k) = (1 / (sigma ^ 2)) * (- lambda(1) * k * ts * cos(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                               - lambda(1) * k * ts * sin(lambda(3) * k * ts + lambda(4)) * Y2(k)...
                                               - lambda(2) * k * ts * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               - lambda(2) * k * ts * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));    
        dd_lambda_53(k) = (1 / (sigma ^ 2)) * (- lambda(2) * k * ts * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               - lambda(2) * k * ts * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));
        
        % Четвертый столбец матрицы
        dd_lambda_14(k) = dd_lambda_41(k);                        
        dd_lambda_24(k) = dd_lambda_42(k);                      
        dd_lambda_34(k) = dd_lambda_43(k);
        dd_lambda_44(k) = (1 / (sigma ^ 2)) * (- lambda(1) * cos(lambda(3) * k * ts + lambda(4)) * Y1(k)...
                                               - lambda(1) * sin(lambda(3) * k * ts + lambda(4)) * Y2(k)...
                                               - lambda(2) * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               - lambda(2) * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));        
        dd_lambda_54(k) = (1 / (sigma ^ 2)) * (- lambda(2) * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               - lambda(2) * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k)); 
        % Пятый столбец матрицы
        dd_lambda_15(k) = dd_lambda_51(k);                        
        dd_lambda_25(k) = dd_lambda_52(k);                      
        dd_lambda_35(k) = dd_lambda_53(k);
        dd_lambda_45(k) = dd_lambda_54(k);
        dd_lambda_55(k) = (1 / (sigma ^ 2)) * (- lambda(2) * cos(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y3(k)...
                                               - lambda(2) * sin(lambda(3) * k * ts + lambda(4) + lambda(5)) * Y4(k));  
        end
    
    d_lambda = [(sum(d_lambda_1) - (M * lambda(1) / (sigma ^ 2))) (sum(d_lambda_2) - (M * lambda(2) / (sigma ^ 2))) sum(d_lambda_3) sum(d_lambda_4) sum(d_lambda_5)];
    dd_lambda = [sum(dd_lambda_11) sum(dd_lambda_12) sum(dd_lambda_13) sum(dd_lambda_14) sum(dd_lambda_15);...
                 sum(dd_lambda_21) sum(dd_lambda_22) sum(dd_lambda_23) sum(dd_lambda_24) sum(dd_lambda_25);...
                 sum(dd_lambda_31) sum(dd_lambda_32) sum(dd_lambda_33) sum(dd_lambda_34) sum(dd_lambda_35);...
                 sum(dd_lambda_41) sum(dd_lambda_42) sum(dd_lambda_43) sum(dd_lambda_44) sum(dd_lambda_45);...
                 sum(dd_lambda_51) sum(dd_lambda_52) sum(dd_lambda_53) sum(dd_lambda_54) sum(dd_lambda_55)];
    
    lambda_old = lambda;
    u_d = d_lambda * inv(dd_lambda);
    lambda = lambda - u_d;
    end
    
    %% Матрица Рыбака
        J44 = M * ((lambda(1) ^ 2) + (lambda(2) ^ 2));
        J55 = M * ((lambda(1) ^ 2) + (lambda(2) ^ 2));

    for k = 1:M
        J33(k) = ((k * ts) ^ 2) * ((lambda(1) ^ 2) + (lambda(2) ^ 2));
        J43(k) = (k * ts) * ((lambda(1) ^ 2) + (lambda(2) ^ 2));
        J53(k) = (k * ts) * ((lambda(2) ^ 2));
    end
    
    J = (1 / (sigma ^ 2)) * [M 0 0 0 0;...
                             0 M 0 0 0;...
                             0 0 sum(J33) sum(J43) sum(J53);...
                             0 0 sum(J43) J44 J55;...
                             0 0 sum(J53) J55 J55];
    %% Граница Крамера-Рао
    D = -inv(J);   
    D_delta_phi = D(5,5)
    
    for k = 1:M
        S1(k) = lambda(1) * cos( lambda(3) * k * ts + lambda(4) );
        S2(k) = lambda(1) * sin( lambda(3) * k * ts + lambda(4) );
        S3(k) = lambda(2) * cos( lambda(3) * k * ts + lambda(4) + lambda(5) );
        S4(k) = lambda(2) * sin( lambda(3) * k * ts + lambda(4) + lambda(5) );
    end  
    
    t = t * 1e6;
    %% Графики 
    figure
    hold on; grid on; grid minor;
    plot(t,Y1,'Color','black','LineWidth',1.5);
    plot(t,Y2,'Color','red','LineWidth',1.5);
    plot(t,Y3,'Color','blue','LineWidth',1.5);
    plot(t,Y4,'Color','green','LineWidth',1.5);
    plot(t,S1,':','Color','black','LineWidth',1.5);
    plot(t,S2,':','Color','red','LineWidth',1.5);
    plot(t,S3,':','Color','blue','LineWidth',1.5);
    plot(t,S4,':','Color','green','LineWidth',1.5);
    %title('title','FontSize',24);
    xlabel('Время, мкс','FontSize',26);
    ylabel('Амплитуда','FontSize',26);
    pause(0.1);
end