function [plt]=plot_orbit(a0,e0,i0,RAAN0,w0,TA0,planet_consts)
% Plots a single orbit around a planet, with eccentricity vector direction in red,
% starting from keplerian elements 
% 
% Inputs:
% a0             : semi-major axis                       (km)
% e0             : eccentricity
% i0             : inclination                           (deg)
% RAAN0          : right ascension of ascending node     (deg)
% w0             : argument of perigee                   (deg)
% TA0            : true anomaly                          (deg)
% planet_consts  : stuct with planet constants
%
% Plots:
% 'Nominal Orbit'
% 
% Functions required:
% - astroConstants()
% - kep2car_rad()
%
% Contributors:
% Gaballo Paolo
% 
% Versions:
% 2021-04-13, third version
% 2021-04-12, second version (plot_orbit_Titan)
% 2021-02-13, first version (plot orbit)

% Constants
% PHYSICAL CONSTANTS
body_mass = planet_consts.mass;                 % Mass of the parent body
Grav_const = 6.67408e-11;                       % Standard Gravitational Constant
R_B = planet_consts.radius;                     % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass;                    % Gravitational Parameter of Body (3.986e14 for earth)

orbit.const.deg     = pi/180;                   %Degrees to radians
orbit.const.mu      = Mu_B;                     %Gravitational parameter (km^3/s^2)
orbit.const.R       = R_B;                      %Planet radius (km)
orbit.const.J2      = astroConstants(9);        %Earth’s J2 %%%%%%%%% CHANGE %%%%%%%%%%%%%
orbit.const.AU      = astroConstants(2);        %Astronomical Unit (AU)

% Conversion from deg to rad
deg     = orbit.const.deg;
i0      = i0 * deg;         %[rad]
RAAN0   = RAAN0 * deg;      %[rad]
w0      = w0 * deg;         %[rad]
TA0     = TA0 * deg;        %[rad]

% Evaluate r0 and v0
[r0, v0] = kep2car_rad(a0,e0,i0,RAAN0,w0,TA0,orbit.const.mu);

% Evaluate orbit period T0 and define number of time steps for ode113
T0      = 2*pi/sqrt(orbit.const.mu)*a0^1.5; %[sec]
t0      = 0;
tf      = T0;
nout    = 4000;
tspan   = linspace(t0, tf, nout);

options = odeset('reltol', 1.e-13, 'abstol', 1.e-14, 'initialstep', T0/1000) ;

[~,xx_car] = ode113(@eq_motion, tspan, [r0 v0], options, orbit);

% Position vector and velocity vector
r_vect = xx_car(:,1:3);
v_vect = xx_car(:,4:6);

% Angular momentum vector
h_vect=cross(r_vect,v_vect);
% Eccentricity vector
e_vect = cross(v_vect,h_vect)/orbit.const.mu - r_vect./sqrt((r_vect(:,1)).^2+(r_vect(:,2)).^2+(r_vect(:,3)).^2);
% Perigee radius
rp = a0*(1-e0);

hold on
% Plot orbit
plt = plot3(r_vect(:,1),r_vect(:,2),r_vect(:,3),'Linewidth',2);
% Plot eccentricity vector
plot3([0,e_vect(1,1)*1.5*rp/e0],[0,e_vect(1,2)*1.5*rp/e0],[0,e_vect(1,3)*1.5*rp/e0],'r','Linewidth',1.5)
xlabel('[km]','FontSize',9)
ylabel('[km]','FontSize',9)
zlabel('[km]','FontSize',9)
title('Nominal Orbit')

end