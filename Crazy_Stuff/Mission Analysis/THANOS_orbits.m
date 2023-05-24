close all
clear all
clc

Titan.mass = 1.3452e23;                          % Mass of the parent body
Grav_const = 6.67408e-11;                       % Standard Gravitational Constant
% c = physconst('LightSpeed');                    % Wave propagation velocity [m/s]
Titan.radius = 2.57473e3;                                % Radius of Body (Earth = physconst('EarthRadius'))
h = 900; %[km]
a = Titan.radius + h;

Titan_3D
p1 = plot_orbit(a,0,45,0,0,0,Titan);
p2 = plot_orbit(a,0,90,0,0,0,Titan);
p3 = plot_orbit(a,0,135,0,0,0,Titan);
legend([p1,p2,p3],'Phase 1 (incl = 45°)','Phase 2 (incl = 90°)','Phase 3 (incl = 135°)')
% D_sat = 1221.865e3; %[km] mean distance of Titan from Saturn
% Sat_mass = 5.6834e26; %[kg]
% Sat_rad = 54364; %[km]
% R_soi = D_sat * (Titan.mass/Sat_mass)^(2/5);
% 
% figure
% Titan_3D
% plot_orbit(5000,0,45,0,0,0,Titan)
% plot_orbit(5000,0,90,0,0,0,Titan)
% plot_orbit(5000,0,135,0,0,0,Titan)
% planet_3D(Sat_rad,D_sat,'saturn_map.jpg')