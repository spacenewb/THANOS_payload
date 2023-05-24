clc
clear all

body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant
J2 = 3.15e-5; % J2 for Titan (Casssini)
ecc = 0;
a = 900e3 + R_B;

Saturn_Year = 365.25*29.4571; % Saturn Year

delta_Omega_SSO = deg2rad(360/Saturn_Year) % Precession of Titan for SSO

%delta_Omega = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(e^2))^2)*a^2)
A = (((1-(ecc^2))^2)*a^2)*delta_Omega_SSO
B = -(3*pi*J2*R_B^2)
C = (A/B)
D = acos(C)

i = rad2deg(acos(((((1-(ecc^2))^2)*a^2)*delta_Omega_SSO)/(-3*pi*J2*R_B^2)));


%%
i = linspace(0,2*pi,1000);
alt_inc = 1000;
a = 900e3 + R_B;
ecc = 0;
delta_Omega1 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); a=a+alt_inc;
delta_Omega2 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); a=a+alt_inc;
delta_Omega3 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); a=a+alt_inc;
delta_Omega4 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); a=a+alt_inc;
delta_Omega5 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); a=a+alt_inc;

figure
plot(rad2deg(i),delta_Omega1,rad2deg(i),delta_Omega2,rad2deg(i),delta_Omega3,rad2deg(i),delta_Omega4,rad2deg(i),delta_Omega5)

%%
i = linspace(0,2*pi,1000);
ecc_inc = 0.15;
a = 900e3 + R_B;
ecc = 0;
delta_Omega1 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); ecc=ecc+ecc_inc;
delta_Omega2 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); ecc=ecc+ecc_inc;
delta_Omega3 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); ecc=ecc+ecc_inc;
delta_Omega4 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); ecc=ecc+ecc_inc;
delta_Omega5 = (-3*pi*J2*cos(deg2rad(i))*R_B^2)/(((1-(ecc^2))^2)*a^2); ecc=ecc+ecc_inc;

figure
plot(rad2deg(i),delta_Omega1,rad2deg(i),delta_Omega2,rad2deg(i),delta_Omega3,rad2deg(i),delta_Omega4,rad2deg(i),delta_Omega5)