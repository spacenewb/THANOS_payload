clear all
close all
clc

%% Orbit Parameters
R_T = 2.57473e3; %[km] Radius of Titan 
h = 900; %[km] orbit altitude

a0 = R_T + h; %[km] semimajor axis
e0 = 0; %[] eccentrity
i0 = 45; %[deg] inclination
RAAN0 = 0; %[deg] right ascension of the ascending node
w0 = 0; %[deg] argument of perigee
TA0 = 0; %[deg] true anomaly
body_mass   = 1.3452e23;                % Mass of the parent body
Grav_const  = 6.67408e-20;              % Standard Gravitational Constant
R_B         = 2.57473e3;                % Radius of Body (Earth = physconst('EarthRadius'))
mu_B        = Grav_const*body_mass;     % Gravitational Parameter of Body 
T0      = 2*pi/sqrt(mu_B)*a0^1.5; %[sec]
omegaT = 2*pi/(15.945421*24*3600);   %[rad/s] 
dlambda = T0*omegaT;
duty_cycle = 0.2;
% n_orbit1 = 2*pi*7/dlambda/duty_cycle
% n_orbit2 = 2*pi/dlambda
n_phase = 3;
swath = 167.164346;
orb_tx = 6;
swath_angle = 2*atan(swath/R_B/2)
n_sw_dlambda =(dlambda/swath_angle)

num_orbit = 2*pi/dlambda

num_phase_orbit = ((orb_tx) / duty_cycle)*( 2*pi/swath_angle)
mission_time = n_phase*num_phase_orbit*T0/(365.25*86400)


% num_phase_orbit_2 = (mission_time *(365.25*86400)/n_phase/T0)
% swath_angle_pi2 = (orb_tx / duty_cycle)/(num_phase_orbit_2)
% swath_angle_2 = swath_angle_pi2*(2*pi)
% swath_2 = round(tan(swath_angle_2/2)*2*R_B,3)
% 
% err = swath_angle - swath_angle_2
%%

Titan_img = imread('TitanTexture.jpg');
pixlen = length(Titan_img);
pix_km = pixlen/(2*pi*R_B);
pix_swath = pix_km * swath


%%

[r_vect,orbit_par] = plot_orbit_Titan(a0,e0,i0,RAAN0,w0,TA0,1.001);%num_phase_orbit/5);
orbit_par = tracks_fun_Titan(r_vect,orbit_par.tspan,orbit_par);
[orbit_par] = plot_GT_Titan(orbit_par,pix_swath,1,0.0005);

% [r_vect,orbit_par] = plot_orbit_Titan(a0,e0,i0,RAAN0,w0,TA0,10);
% orbit_par = tracks_fun_Titan(r_vect,orbit_par.tspan,orbit_par);
% [orbit_par] = plot_GT_Titan(orbit_par,1,1,1);
%%
% R_T = 2.57473e3; %[km] Radius of Titan 
% h = 900; %[km] orbit altitude
% 
% a0 = R_T + h; %[km] semimajor axis
% e0 = 0; %[] eccentrity
% i0 = 89.5; %[deg] inclination
% RAAN0 = 0; %[deg] right ascension of the ascending node
% w0 = 0; %[deg] argument of perigee
% TA0 = 0; %[deg] true anomaly
% n_orbit=100;
% 
% 
% [r_vect,orbit_par] = plot_orbit_Titan(a0,e0,i0,RAAN0,w0,TA0,n_orbit);
% orbit_par = tracks_fun_Titan(r_vect,orbit_par.tspan,orbit_par);
% [orbit_par] = plot_GT_Titan(orbit_par);