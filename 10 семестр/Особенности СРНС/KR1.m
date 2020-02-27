close all; clear all; clc;

Time_year = 2014;
Time_month = 2;
Time_day = 28;

Time_hour = 12;
Time_minutes = 0;
Time_seconds = 0;


N4 = floor(1+(Time_year-1996)/4);
N_T = 365*(Time_year-1996-4*(N4-1)) + 31 + Time_day;
t = Time_seconds + Time_minutes*60 + Time_hour*60*60 + 10800;

Time_GLN = [N4 N_T t];

%t_i = 
%t_labda_A =
%N = 
%N_A = 
%delta_t_pr = delta_N_A*86400 + (t_i-t_labda_A)