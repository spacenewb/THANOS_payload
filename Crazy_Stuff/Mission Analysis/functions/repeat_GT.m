function [a,m,orbit] = repeat_GT(orbit,k,pert)
% Computes the semimajor axis for a repeating orbit for unperturbed case
% or for J2 perturbed case using fsolve 
%
% Inputs:
% orbit     : structure with data and constants
% k         : number of spacecraft revolutions
% m         : number of Earth rotations
% pert      : J2 perturbed case or not (true or false)
%
% Outputs:
% orbit     : input structure with new data (repeating semi-major axis)
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version 

% k*dlamba = m*2pi
% k = m*pi/dlamda;
% n = 2*pi/T0;
% we/n = dlamba/pi
T0 = orbit.const.T0;
omegaT              = orbit.omegaT;
dlamba = T0*omegaT
m = k*dlamba/2/pi



%omegaE              = deg2rad(15.04)/3600;          %Earth's angular velocity [rad]

mu                  = orbit.const.mu;
R                   = orbit.const.R;
J2                  = orbit.const.J2;
orbit.repeat.k      = k;                            %s/c revolutions
orbit.repeat.m      = m;                            %Earth rotations
a                   = orbit.kep0.a;                 %Semimajor axis (km)
e                   = orbit.kep0.e;                 %eccentricity
RAAN                = rad2deg(orbit.kep0.RAAN);     %Right ascencion of the node (radians)
i                   = rad2deg(orbit.kep0.i);        %Inclination (radians)
w                   = rad2deg(orbit.kep0.w);        %Argument of perigee (radians)
TA                  = rad2deg(orbit.kep0.TA);       %True anomaly (radians)

if pert == false
    %without consider perturbation
    orbit.repeat.a = (mu*m^2/(k*omegaT)^2)^(1/3);   %New semimajor axis for repeated orbit
    a = orbit.repeat.a
elseif pert == true
    %with J2 perturbation effects
    n0 = sqrt(mu/a^3);                                  %mean anomaly
    
    % W_dot0 = -(3/2*(sqrt(mu)*J2*Re^2)/(((1-ecc^2)^2)*a^(7/2)))*cosd(i);
    % w_dot0 = -(3/2*(sqrt(mu)*J2*Re^2)/(((1-ecc^2)^2)*a^(7/2)))*((5/2)*(sind(i))^2-2);
    % M_dot0 = -(3/2*(sqrt(mu)*J2*Re^2)/(((1-ecc^2)^2)*a^(7/2)))*(1-(3/2)*(sind(i))^2);
    
    fun = @(g) [ g(2)+(3/2*(sqrt(mu)*J2*R^2)/(((1-e^2)^2)*g(1)^(7/2)))*cosd(i);
        g(3)+(3/2*(sqrt(mu)*J2*R^2)/(((1-e^2)^2)*g(1)^(7/2)))*((5/2)*(sind(i))^2-2);
        g(4)+(3/2*(sqrt(mu)*J2*R^2)/(((1-e^2)^2)*g(1)^(7/2)))*(1-(3/2)*(sind(i))^2);
        g(5)-sqrt(mu/g(1)^3);
        (omegaE-g(2))/(g(5)+g(3)+g(4))-m/k];
    
    sol0 = [a,RAAN,w,TA,n0];     
    sol = fsolve(fun,sol0);      
    
    orbit.repeat.a = sol(1);
end
end