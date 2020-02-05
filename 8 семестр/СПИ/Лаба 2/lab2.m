close all; 
clear all;
%%
CKO = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2];
GLAZ1 = [1000 977 938 837 845 809 735 667 625 510 365 490 453 239 382 327 330 214 280 82 0]/1000;
N = 2;
coeff1 = polyfit(CKO, GLAZ1, N);
x1 = min(CKO):0.1:max(CKO);
y1 = 0;
for k=0:N
    y1 = y1 + coeff1(N-k+1) * x1.^k;
end

figure;
plot (CKO, GLAZ1,'*','Color','black','LineWidth',1.5);
hold on;
plot (x1, y1,'Color','black','LineWidth',1.5);
grid on;
grid minor;
%title('Зависимость величины раскрыва глаза от СКО шума','FontSize',24);
xlabel('СКО шума, В','FontSize',26);
ylabel('Величина раскрыва глаза, В','FontSize',26);
%%
DeltaFi = [0 15 30 45 60 75 90];
GLAZ2 = [1000 964 868 710 500 260 0]/1000;
N = 2;
coeff2 = polyfit(DeltaFi, GLAZ2, N);
x2 = min(DeltaFi):15:max(DeltaFi);
y2 = 0;
for k=0:N
    y2 = y2 + coeff2(N-k+1) * x2.^k;
end

figure;
plot (DeltaFi, GLAZ2,'*','Color','black','LineWidth',1.5);
hold on;
plot (x2, y2,'Color','black','LineWidth',1.5);
grid on;
grid minor;
%title('Зависимость величины раскрыва глаза от расфазирования','FontSize',24);
xlabel('Расфазирование, град','FontSize',26);
ylabel('Величина раскрыва глаза, В','FontSize',26);
%%
CKO = [1 0.85 0.66 0.42 0];
DeltaFi = [0 15 30 45 60];
N = 2;
coeff3 = polyfit(DeltaFi, CKO, N);
x3 = min(DeltaFi):15:max(DeltaFi);
y3 = 0;
for k=0:N
    y3 = y3 + coeff3(N-k+1) * x3.^k;
end

figure;
plot (DeltaFi,CKO,'*','Color','black','LineWidth',1.5);
hold on;
plot (x3, y3,'Color','black','LineWidth',1.5);
grid on;
grid minor;
%title('Диаграмма обмена','FontSize',24);
xlabel('Расфазирование, град','FontSize',26);
ylabel('СКО шума, В','FontSize',26);

%%
P = [0 3.8e-4 1.1e-2 0.7e-1 2.4e-1 3.2e-1];
DeltaFi = [0 2 3 5 10 15];
N = 5;
coeff4 = polyfit(DeltaFi, P, N);
x4 = min(DeltaFi):1:max(DeltaFi);
y4 = 0;
for k=0:N
    y4 = y4 + coeff4(N-k+1) * x4.^k;
end

figure;
semilogy(DeltaFi,P,'*','Color','black','LineWidth',1.5);
hold on;
semilogy(x4, y4,'Color','black','LineWidth',1.5);
grid on;
grid minor;
xlim ([0 15]);
ylim ([1e-4 1]);
%title('Диаграмма обмена','FontSize',24);
xlabel('Расфазирование, град','FontSize',26);
ylabel('Вероятность ошибки','FontSize',26);

