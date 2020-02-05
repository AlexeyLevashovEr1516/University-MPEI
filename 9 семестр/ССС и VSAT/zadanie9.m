function Main()
    close all; clear all; clc;
    
    M = 1;
    Rvh = 3*1e6;
    f = 1*1e6;
    
    t = 0:1e-9:8/f;
    
    for k = 1:length(t)
        kt = k*1e-9;
        if kt < 1*M/f
            a(k) = 0/7;
            code(k)
        elseif (kt >= 1*M/f) && (kt < 2*M/f)
            a(k) = 1/7;
        elseif (kt >= 2*M/f) && (kt < 3*M/f)
            a(k) = 2/7;
        elseif (kt >= 3*M/f) && (kt < 4*M/f)
            a(k) = 3/7;
        elseif (kt >= 4*M/f) && (kt < 5*M/f)
            a(k) = 4/7;
        elseif (kt >= 5*M/f) && (kt < 6*M/f)
            a(k) = 5/7;
        elseif (kt >= 6*M/f) && (kt < 7*M/f)
            a(k) = 6/7;
        elseif kt >= 7*M/f
            a(k) = 7/7;
        end
        y(k) = a(k)*sin(4*pi*f*kt);
    end
    
    figure(1);
%     subplot(3,1,1)
%     plot(t*1e6, bit, 'LineWidth', 2)
%     hold on, grid on;
%     ylim([0, 1])
%     xlabel("Время, мкс");
%     ylabel("Двоичная последовательность");
%     
%     subplot(3,1,2)
%     plot(t*1e6, phi*180/pi, 'LineWidth', 2)
%     hold on, grid on;
%     ylim([-180, 180])
%     xlabel("Время, мкс");
%     ylabel("Начальная фаза, град");
%     
%     subplot(3,1,3)
    plot(t*1e6, y, 'LineWidth', 2)
    hold on, grid on;
    plot(t*1e6, [a; -a], ':', 'LineWidth', 2)
    xlabel("Время, мкс");
    ylabel("Амплитуда, В");

end