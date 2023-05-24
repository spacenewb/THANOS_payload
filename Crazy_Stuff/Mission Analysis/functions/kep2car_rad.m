function [r, v] = kep2car_rad(a,ecc,i,W,w,theta, mu)
% Returns position vector and velocity vector starting from keplerian elements in input
% All the angles must be in rad
%
% Example:
% [r, v] = kep2car_rad(a,ecc,i,W,w,theta) : mu is setted as default
%
% Inputs:
% a         : semi-major axis                       (km)
% ecc       : eccentricity
% i         : inclination                           (deg)
% W         : right ascension of ascending node     (deg)
% w         : argument of perigee                   (deg)
% theta     : true anomaly                          (deg)
% mu        : gravitational parameter               (km^3/s^2)
%
% Outputs:
% r         : position vector                       (km)   (size(1,3))
% v         : velocity vector                       (km/s) (size(1,3))
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version

if nargin == 6
    mu = 398600.433;
end
%matrici di rotazione;
R1 = [cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3 = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];
T = R1'*R2'*R3';

p = a*(1-ecc^2);
r = p /(1+ ecc*cos(theta));
r = r*[cos(theta) sin(theta) 0]';
v = sqrt(mu/p)*[-sin(theta) ecc+cos(theta) 0]';

r = T*r;
v = T*v;
