close all
clear all
clc
% Data volume (Housekeeping - Orbital - Attitude) minimum case
Titan_rate = 4e6; %430e3; %[bit/sec]

%--------S/C Orbit--------------%
R_t=2.57473e3; %[Km]
% mu_t=8971.15; %[Km^3/s^2]
body_mass   = 1.3452e23;                % Mass of the parent body
Grav_const  = 6.67408e-20;              % Standard Gravitational Constant
R_B         = 2.57473e3;                % Radius of Body (Earth = physconst('EarthRadius'))
mu_T        = Grav_const*body_mass;     % Gravitational Parameter of Body 
semajax_sc=900+R_t; %[Km]
ecc_sc=0; %[/]
period_sc= 2*pi*sqrt(semajax_sc^3/mu_T); %[s]
incl_sc=deg2rad(45); %[]
mean_motion_sc= 2*pi/period_sc;
% T_orbit = 5.5; %[h]
T_orbit = period_sc/3600; %[h]

%% Housekeeping data
% Housekeeping sampling frequency is very low (T=1 sec - 2 min),
% consequently:
min_bwidth = 10; %[bit/s]
typ_bwidth = 100; %[bit/s]
max_bwidth = 10*1024*8; %[bit/s]

% Temperature: (P/L, batteries, solar arrays, Thrusters, EPS)
range_Temp = 400 - 0; %[K]
n_Temp=1;
while 2^n_Temp < range_Temp
    n_Temp = n_Temp+1;
end
%single data size (48 bits for Source Packet Header + 16 bits for Data
%Field Header (ancillary data, s/c time))
packet_Temp = 48+16+n_Temp+1; %[bit]
vol_Temp = packet_Temp*T_orbit*60; %[bit]
% 5 differents compponents to monitor
vol_Temp_est = vol_Temp*5; %[bit]
vol_Temp_typ = typ_bwidth * T_orbit * 3600; %[bit]
vol_Temp = max(vol_Temp_est, vol_Temp_typ);

% Pressure, structural stress, flow rate: (Fuel tanks, structures)
vol_pres = typ_bwidth * T_orbit * 3600; %[bit]

% Voltages and currents:
vol_volt = typ_bwidth * T_orbit * 3600; %[bit]

% Operative modes: (P/L status, Heaters) (we can consider min bandwidth)
vol_modes = typ_bwidth * T_orbit * 3600; %[bit]

% Actuation/deployment of mechanisms (we can consider min bandwidth)
vol_mech = typ_bwidth * T_orbit * 3600; %[bit]

% Total data volume for Housekeeping
Tot_House = vol_Temp + vol_pres + vol_volt + vol_modes + vol_mech ; %vol_health %[bit]

Tot_House_Mbits = Tot_House/(1024^2) %[MBytes]

% House_rate = Tot_House/(T_orbit*3600); %[bit/sec]

%% Payload Health check data:
vol_health = typ_bwidth * T_orbit * 3600; %[bit]

Tot_Health_Mbits = vol_health/(1024^2) %[MBytes]

%% Orbital data
% Position s/c (each second)
range_pos = 50000 - 0; %[cm]
n_pos=1;
while 2^n_pos < range_pos
    n_pos = n_pos+1;
end
%single data size (48 bits for Source Packet Header + 16 bits for Data
%Field Header (ancillary data, s/c time))
packet_pos = 48+16+n_pos+1; %[bit]
vol_pos = 3*packet_pos*T_orbit*3600; %[bit] (3 components)

% Velocity s/c (each second)
sc_vel = 1500; %[m/s] about
range_vel = 40000 - 0; %[cm/s]
n_vel=1;
while 2^n_vel < range_vel
    n_vel = n_vel+1;
end
%single data size (48 bits for Source Packet Header + 16 bits for Data
%Field Header (ancillary data, s/c time))
packet_vel = 48+16+n_vel+1; %[bit]
vol_vel = 3*packet_vel*T_orbit*3600; %[bit] (3 components)

Tot_Orbital = vol_pos + vol_vel; %[bit]
Tot_Orbital_Mbits = Tot_Orbital/(1024^2) %[MBytes]

% Orbital_rate = Tot_Orbital/(T_orbit*3600); %[bit/sec]

%% Attitude data
% generaldata rate in the range of 10-1k Bytes/sec
% consider the worst case
ACS_bwidth = 8*1024; %[bit/s]
vol_att = ACS_bwidth * T_orbit * 3600; %[bit]

Tot_Attitude = vol_att; %[bit]
Tot_Attitude_Mbits = Tot_Attitude/(1024^2) %[MBytes]

% Attitude_rate = Tot_Attitude/(T_orbit*3600); %[bit/sec]

%% Total

TOT_vol_typical = Tot_House + Tot_Orbital + Tot_Attitude + vol_health; %[bit]
TOT_vol_typical_Mbits = TOT_vol_typical/(1024^2) %[MBytes]

TOT_data_rate = TOT_vol_typical/(T_orbit*3600); %[bit/sec]

%Transmission Time needed
Trans_time_ = TOT_vol_typical/Titan_rate %[sec] (Titan_rate = 4Mbps)

Data_per_10sec = TOT_vol_typical_Mbits/(T_orbit*360)
Data_per_20sec = TOT_vol_typical_Mbits/(T_orbit*180)
