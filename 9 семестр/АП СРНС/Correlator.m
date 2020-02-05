clear all; clc; close all;

N = 30;

Fd = 99.375e6;
Td = 1/Fd;
T = 0.001;
fif = 8.34e6;

fd     = 100 * N;
fd_rep = fd;

phi     = 10 * N;
phi_rep = 20 * N;

qcno_dB = 27 + N;
qcno    = 10^(qcno_dB/10);

stdn = 50/3;
A    = sqrt(qcno * Td) * 2 * stdn; % qcno = A^2 / (4 stdn^2 Td)

fprintf('A = %f, stdn = %f\n', A, stdn);

L = fix(T * Fd);
t = (0:(L-1)) * Td;

PRN_Length = 1023;
PRN = (-1).^(rand(1, PRN_Length) > 0.5);

tau = 100500e-6;
tau_rep = tau;

nchip     = mod(fix(PRN_Length*(t - tau    )/T), PRN_Length) + 1;
nchip_rep = mod(fix(PRN_Length*(t - tau_rep)/T), PRN_Length) + 1;

Gc     = PRN(nchip    );
Gc_rep = PRN(nchip_rep);

So        = A * Gc     .* cos(2*pi*fif*t + 2*pi*fd    *t + deg2rad(phi)    );
S_rep_cos =     Gc_rep .* cos(2*pi*fif*t + 2*pi*fd_rep*t + deg2rad(phi_rep));
S_rep_sin =     Gc_rep .* sin(2*pi*fif*t + 2*pi*fd_rep*t + deg2rad(phi_rep));

K = 1000; I = nan(1, K); Q = nan(1,K);
for k = 1:K
    n = randn(1,L)*stdn;
  
    Gd     = (-1).^(rand(1,1) > 0.5);
    S = Gd * So;
    
    y = S + n;
    
    I(k) = y * S_rep_cos';
    Q(k) = y * S_rep_sin';
end

maxIQ = 1.1*max(abs(I + 1i*Q));

figure(1);
plot(t*1e3, [y; S; A*Gc]);
xlabel('t, ms')
ylabel('y, S');
grid on

figure(2);
plot(I, Q, '*')
hold on
plot([-maxIQ maxIQ], [0 0], 'k');
plot([0 0], [-maxIQ maxIQ], 'k');
quiver(-maxIQ, 0, 2*maxIQ, 0, 1, 'k');
quiver(0, -maxIQ, 0, 2*maxIQ, 1, 'k');
hold off
xlabel('I')
ylabel('Q');
axis equal
grid on