clc
clear all

% body_mass = 1.3452e23; % Mass of the parent body
% Grav_const = 6.67408e-11; % Standard Gravitational Constant
% R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
% Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
% 
% J2 = 3.15e-5; % J2 for Titan (Casssini)
% 
% ecc = 0;
% a = 900e3 + R_B;
% inc = 120;
% 
% Saturn_Year = 365.25*29.4571; % Saturn Year
% Target_delta_Omega = deg2rad(360/Saturn_Year) % Precession of Titan for SSO
% 
% delta_Omega = (-3*pi*J2*cos(deg2rad(inc))*R_B^2)/(((1-(ecc^2))^2)*a^2)
R_B = 2.57473e6;

lb = [0, (R_B+900e3)]
ub = [359.99, (R_B+2000e3)]

SunSync(0,R_B+900e3)

opt_inc = ga(@SunSync,2,[],[],[],[],lb,ub)

function[error_delta_Omega] = SunSync(inc,a)
    R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))

    J2 = 3.15e-5; % J2 for Titan (Casssini)

    ecc = 0;

    Saturn_Year = 365.25*29.4571; % Saturn Year
    Target_delta_Omega = deg2rad(360/Saturn_Year); % Precession of Titan for SSO
    
    delta_Omega = (-3*pi*J2*cos(deg2rad(inc))*R_B^2)/(((1-(ecc^2))^2)*a^2);
    
    error_delta_Omega = Target_delta_Omega - delta_Omega;
end