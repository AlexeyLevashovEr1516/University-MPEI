close all; 
clear all;

%%
f = [995 996 997 998 999 1000 1001 1002 1003 1004 1005];
Gain = [10 8 6 4 2 0.2 2 4 6 8 10];

%N = 8;
%coeff1 = polyfit(f, Gain, N);
%x1 = min(f):1:max(f);
%y1 = 0;
%for k=0:N
%    y1 = y1 + coeff1(N-k+1) * x1.^k;
%end


figure;
plot(f,Gain,'*','Color','black','LineWidth',1.5);
hold on;
%plot(x1,y1,'Color','black','LineWidth',1.5);
grid on;
grid minor;
%xlim ([1e-4 1]);
%ylim ([1e-4 1]);
title('Зависимость полосы захвата СВН при ограничении время захвата 10с','FontSize',24);
xlabel('Полоса захвата f, Гц','FontSize',26);
ylabel('Коэффициент усиления Gain','FontSize',26);

