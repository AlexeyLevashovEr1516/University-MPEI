function Main()
    close all; clear all; clc;
    
    T = 10e-3;
    alpha = 1;
    c = 3e8;
    f0 = 1602e6;
    omega0 = 2*pi*f0;
    q_dB = 30;
    q = 10^(q_dB/10);
    N0 = 2/(q*T^2)*(1+1/(2*q*T));
    
    sigma_a = 1:30;
      
    sigma_OMEGA = nan(size(sigma_a));
    
    for k = 1:30        
        S_ksi = 2*k^2*alpha*(omega0/c)^2;
        K1 = alpha*(sqrt(1+2*sqrt(S_ksi/(alpha^2*N0)))-1);
        K2 = (K1^2)/2;
        D11 = K1*N0/2;
        sigma_OMEGA_k = sqrt(D11);
        sigma_OMEGA(k) = sigma_OMEGA_k;
        deltaF_k = (K1^2*K2+(alpha*K1+K2)^2)/(K2*(K1+alpha));
        deltaF(k) = deltaF_k;
    end 
    
    figure(1);
    hold on; grid on;
    plot(sigma_a, sigma_OMEGA);
    xlabel('Sigma a, m/s^2');
    ylabel('Sigma Omega, Hz');
    
    figure(2);
    hold on; grid on;
    plot(sigma_a, deltaF);
    xlabel('Sigma a, m/s^2');
    ylabel('deltaF, Hz');
end