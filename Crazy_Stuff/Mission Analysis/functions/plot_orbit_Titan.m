function [r_vect,orbit] = plot_orbit_Titan(a0,e0,i0,RAAN0,w0,TA0,n_orbit)
% Plots a single orbit around Titan, with eccentricity vector direction in red,
% starting from keplerian elements 
% 
% Inputs:
% a0             : semi-major axis                       (km)
% e0             : eccentricity
% i0             : inclination                           (deg)
% RAAN0          : right ascension of ascending node     (deg)
% w0             : argument of perigee                   (deg)
% TA0            : true anomaly                          (deg)
% n_orbit        : number of orbit   
%
% Outputs:
% r_vect         : position vector                       (km)
% orbit          : structure with all the orbit data     (km)
%
% Plots:
% 'Nominal Orbit'
% 
% Functions required:
% - astroConstants()
% - kep2car_rad()
% - Titan_3D()
%
% Contributors:
% Gaballo Paolo
% 
% Versions:
% 2021-04-13, second version
% 2021-02-13, first version (plot orbit)

% Constants
% PHYSICAL CONSTANTS
body_mass   = 1.3452e23;                % Mass of the parent body
Grav_const  = 6.67408e-20;              % Standard Gravitational Constant
R_B         = 2.57473e3;                % Radius of Body (Earth = physconst('EarthRadius'))
mu_B        = Grav_const*body_mass;     % Gravitational Parameter of Body (3.986e14 for earth)
% omegaT      = 2*pi/(((15*24+22)*60+41)*60)   % Angular velocity of Titan [rad/s]
omegaT      = 2*pi/(15.945421*24*3600);  % Angular velocity of Titan [rad/s]
orbit.const.deg     = pi/180;           % Degrees to radians
orbit.const.mu      = mu_B;             % Gravitational parameter (km^3/s^2)
orbit.const.R       = R_B;              % Titan’s radius (km)
orbit.const.J2      = 33.089e-6;        % Titan’s J2 (36.7e-6 for Cassini)

% Conversion from deg to rad
deg     = orbit.const.deg;
i0      = i0 * deg;         %[rad]
RAAN0   = RAAN0 * deg;      %[rad]
w0      = w0 * deg;         %[rad]
TA0     = TA0 * deg;        %[rad]

orbit.kep0.a = a0;
orbit.kep0.i = i0;
orbit.kep0.e = e0;
orbit.kep0.RAAN = RAAN0;
orbit.kep0.w = w0;
orbit.kep0.TA = TA0;
orbit.omegaT = omegaT;

% Evaluate r0 and v0
[r0, v0] = kep2car_rad(a0,e0,i0,RAAN0,w0,TA0,orbit.const.mu);

% Evaluate orbit period T0 and define number of time steps for ode113
T0      = 2*pi/sqrt(orbit.const.mu)*a0^1.5; %[sec]
orbit.const.T0 = T0;
t0      = 0;
tf      = n_orbit*T0;
nout    = 4000*n_orbit;
tspan   = linspace(t0, tf, nout);
orbit.tspan = tspan;

options = odeset('reltol', 1.e-13, 'abstol', 1.e-14, 'initialstep', T0/1000) ;

[~,xx_car] = ode113(@eq_motion, tspan, [r0 v0], options, orbit,'J2');

% Position vector and velocity vector
r_vect = xx_car(:,1:3);
v_vect = xx_car(:,4:6);

% Angular momentum vector
h_vect = cross(r_vect,v_vect);
% Eccentricity vector
e_vect = cross(v_vect,h_vect)/orbit.const.mu - r_vect./sqrt((r_vect(:,1)).^2+(r_vect(:,2)).^2+(r_vect(:,3)).^2);
% Perigee radius
rp = a0*(1-e0);

figure
% Plot Titan
Titan_3D()
hold on
% Plot orbit
plot3(r_vect(:,1),r_vect(:,2),r_vect(:,3),'Linewidth',1);
% Plot eccentricity vector
plot3([0,e_vect(1,1)*1.5*rp/e0],[0,e_vect(1,2)*1.5*rp/e0],[0,e_vect(1,3)*1.5*rp/e0],'r','Linewidth',1.5)
xlabel('[km]','FontSize',9)
ylabel('[km]','FontSize',9)
zlabel('[km]','FontSize',9)
title('Nominal Orbit')

end