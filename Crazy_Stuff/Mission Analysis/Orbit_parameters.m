clear all
close all
clc

% WS= 14.82; %[W/m^2] Solar irradiance at Saturn
% Tss= -138.5; %[deg C] Saturn surface temperature at 1 bar
% WIR=18.2828; %[W/m^2]
% BBT= 81; %[K] Black Body Temperature
% Bond_alb= 0.342; %[/] Bond Albedo
% Geom_alb= 0.499; %[/] Geometrica Albedo
% Sat_Surf= 4.2*10^16; %[m^2] Surface of Saturn
% Tds= 2.7; %[K] Deep Space Temperauture
diam_s=120536; %[Km]
R_s=diam_s/2;  %[Km]
R_t=2576; %[Km]
t_simul=2*365*24*3600; %[s]


%---------Saturn Orbit-----------%
mu_sun=132712440018; %[Km^3/s^2]
peri_s=1352550000;    %[Km]
apo_s=1515500000;      %[Km]
semajax_s=(apo_s+peri_s)/2;  %[Km]
ecc_s= (apo_s-peri_s)/(apo_s + peri_s); %[/]
% period_s=29.45*365*24*3600; %[s]
period_s = 2*pi*sqrt(semajax_s^3/mu_sun);
mean_motion_s= (2*pi)/period_s; %[rad/s]
incl_s=deg2rad(2.485); %[rad]
incl_ax_s=deg2rad(26.73); % [rad]
mean_r_s=1.427*10^9; %[Km]
% t_simul=period_s/2;

%---------Earth Orbit-----------%
% mu_sun=132712440018; %[Km^3/s^2]
% peri_s=1352550000;    %[Km]
% apo_s=1515500000;      %[Km]
semajax_E= 149.596e6; %(apo_s+peri_s)/2;  %[Km]
ecc_E= 0.01671022; %(apo_s-peri_s)/(apo_s + peri_s); %[/]
% period_s=29.45*365*24*3600; %[s]
period_E = 2*pi*sqrt(semajax_E^3/mu_sun);
mean_motion_E= (2*pi)/period_s; %[rad/s]
incl_E=deg2rad(0.00005); %[rad]
incl_ax_E=deg2rad(23.44); % [rad]
mean_r_E=semajax_s; %[Km]
% t_simul=period_s/2;

%--------Titan Orbit------------%
mu_s=37931187; %[Km^3/s^2]
semajax_t=1221830; %[Km]
ecc_t=0.0293; %[/]
% period_t=15.945421*24*3600; %[s]
period_t = 2*pi*sqrt(semajax_t^3/mu_s);
incl_t=deg2rad(0.3485); %[rad]
mean_motion_t= (2*pi)/period_t; %[rad/s]
incl_t_ecl=incl_t+incl_s;

%--------S/C Orbit--------------%
mu_t=8971.15; %[Km^3/s^2]
semajax_sc=1500+R_t; %[Km]
ecc_sc=0; %[/]
period_sc= 2*pi*sqrt(semajax_sc^3/mu_t); %[s]
incl_sc=deg2rad(45); %[]
mean_motion_sc= 2*pi/period_sc;

%max_step is the integration max step to input in simulink when the 
%integrator is characterized. It is important to set this value small 
%enough wrt the integration time otherwhise the eclipse block could miss 
%some points of integration, thus some eclipses are not computed. 
max_step=t_simul/8.5e3;
%------------PROVA--------------%
% M_earth= 5.97219*10^24;
% R_E=6371;
% G=6.6743015*10^-20; 
% mu=(G*M_earth);
% mu_s=132712440018;
% 
% a=15000+R_E;
% ecc=0;
% T_orb=((2*pi))*(((a^3)/(mu))^0.5);
% n_orb=(2*pi)/T_orb;
% m=n_orb;
% 
% 
% incl=deg2rad(0);
% incl_E=deg2rad(0);
% AU=1.495978707*10^8;
% T_E=365*(24*60*60);
