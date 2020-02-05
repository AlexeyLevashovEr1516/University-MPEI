close all; 
clear all;

%%
P = [1e-4 0.5e-3 1e-3 0.5e-2 1e-2 2e-2 3e-2 0.5e-1 1e-1 0.5];
Pdsk = [1.8e-4 1.8e-4 3.6e-4 0.58e-2 1.03e-2 2.2e-2 3e-2 0.54e-1 0.96e-1 0.49];
Pdec = [0 0 0 0 1e-10 2.5e-3 1.2e-2 4.7e-2 1e-1 0.48];

N = 1;
coeff1 = polyfit(P, Pdsk, N);
x1 = min(P):1e-4:max(P);
y1 = 0;
for k=0:N
    y1 = y1 + coeff1(N-k+1) * x1.^k;
end

N = 1;
coeff2 = polyfit(P, Pdec, N);
x2 = min(P):1e-4:max(P);
y2 = 0;
for k=0:N
    y2 = y2 + coeff2(N-k+1) * x2.^k;
end

figure;
loglog(P,Pdsk,'*','Color','black','LineWidth',1.5);
hold on;
loglog(x1,y1,'Color','black','LineWidth',1.5);
hold on;
loglog(P,Pdec,'*','Color','red','LineWidth',1.5);
hold on;
loglog(x2,y2,'Color','red','LineWidth',1.5);
grid on;
grid minor;
xlim ([1e-4 1]);
ylim ([1e-4 1]);
%title('Диаграмма обмена','FontSize',24);
%xlabel('Расфазирование, град','FontSize',26);
%ylabel('Вероятность ошибки','FontSize',26);

