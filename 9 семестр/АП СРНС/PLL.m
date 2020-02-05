clear all; clc; close all;

T = 0.010;

Xist = [10; 
         1];
     
Xfilt = [0;
         0];

Xfilt_FLL = [0;
             0];
     
F = [1 T;
     0 1];

Tmod = 10;
K = fix(Tmod / T);
 
G = [0;
     1];

std_ksi = 1.3 * T;
ksi = std_ksi * randn(1, K);

qcno_dB = 45;
qcno = 10^(qcno_dB/10);

stdnIQ = 7;
nI = stdnIQ * randn(1, K);
nQ = stdnIQ * randn(1, K);
A_IQ = sqrt(2*qcno*T) * stdnIQ;  

dF_PLL = 20;
Kfilt_PLL = [ 8/3 * dF_PLL   * T;
              32/9 * dF_PLL^2 * T];

dF_FLL = 3;
Kfilt_FLL = [ 8/3 * dF_FLL   * T;
              32/9 * dF_FLL^2 * T];
          
t = (1:K)*T;
phi_ist  = nan(1, K);
  w_ist  = nan(1, K);
phi_filt = nan(1, K);
  w_filt = nan(1, K);
phi_filt_FLL = nan(1, K);
  w_filt_FLL = nan(1, K);
  
phi_corr = 0; Iold = 0; Qold = 0;
for k = 1:K
    Xist = F * Xist + G * ksi(k);
    
    phi_ist(k) = Xist(1);
      w_ist(k) = Xist(2);
    
    % PLL
    Xfilt_extr = F*Xfilt;
    phi_extr = Xfilt_extr(1);
    dPhi = phi_ist(k) - phi_extr;

    Q = - A_IQ * sin(dPhi) + nQ(k);
    Ud = -Q;  
    Sd = A_IQ;  
    
    Xfilt = Xfilt_extr + Kfilt_PLL * Ud/Sd;

    phi_filt(k) = Xfilt(1);
      w_filt(k) = Xfilt(2);
    
    % FLL
    Xfilt_extr_FLL = F*Xfilt_FLL;
    w_extr = Xfilt_extr_FLL(1);
    dPhi = phi_ist(k) - phi_corr;
    phi_corr = phi_corr + w_extr * T;
    dW = w_ist(k) - w_extr;

    I = - A_IQ * cos(dPhi) + nI(k);
    Q = - A_IQ * sin(dPhi) + nQ(k);
    Ud = Q*Iold - I*Qold;  
    Sd = A_IQ;  
    Iold = I; Qold = Q;
    
    Xfilt_FLL = Xfilt_extr_FLL + Kfilt_FLL * Ud/Sd;

    phi_filt_FLL(k) = phi_corr;
      w_filt_FLL(k) = Xfilt_FLL(1);
    
    
end

figure(1)
subplot(2,1,1)
plot(t, phi_ist/2/pi, t, phi_filt/2/pi, t, phi_filt_FLL/2/pi)
ylabel('\phi, cycles')
subplot(2,1,2)
plot(t, rad2deg(w_ist), t, rad2deg(w_filt), t, rad2deg(w_filt_FLL))
ylabel('\omega, Hz')
xlabel('t, sec')