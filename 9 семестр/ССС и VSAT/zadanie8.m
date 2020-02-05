function Main()
    close all; clear all; clc;
    
    N = 3;
    M = 10/7;
    f = 10/7*1e6;   
    t = 0:1e-9:N*M/f;

    phi0 = 33.75*pi/180;
    phi1 = -95.625*pi/180;
    phi2 = 129.375*pi/180;
    
    code = [1 0 0 1 1 0 0 0 1 1 1 1 1 1 0 1 1 1 0 0 0 0 0 0 0 0];
    num = 1;
    for k = 1:length(t)
        kt = k*1e-9;
        bit(k) = code(num);
        if kt >= M*num/f/6
            num = num + 1;
        end
    end
    
    for k = 1:length(t)
        kt = k*1e-9;
        if kt < 1*M/f
            y(k) = cos(2*pi*f*kt - phi0);
            phi(k) = phi0;
        elseif (kt >= 1*M/f) && (kt < 2*M/f)
            y(k) = cos(2*pi*f*kt - phi1);
            phi(k) = phi1;
        elseif kt >= 2*M/f
            y(k) = cos(2*pi*f*kt - phi2);
            phi(k) = phi2;
        end
    end
    
    figure(1);
    subplot(3,1,1)
    plot(t*1e6, bit, 'LineWidth', 2)
    hold on, grid on;
    ylim([0, 1])
    xlabel("Время, мкс");
    ylabel("Двоичная последовательность");
    
    subplot(3,1,2)
    plot(t*1e6, phi*180/pi, 'LineWidth', 2)
    hold on, grid on;
    ylim([-180, 180])
    xlabel("Время, мкс");
    ylabel("Начальная фаза, град");
    
    subplot(3,1,3)
    plot(t*1e6, y, 'LineWidth', 2)
    hold on, grid on;
    xlabel("Время, мкс");
    ylabel("Амплитуда, В");

end