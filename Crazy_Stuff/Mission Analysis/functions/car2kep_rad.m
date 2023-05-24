function [a,ecc,i,W,w,theta] = car2kep_rad(r,v,mu)
% Takes in input position vector r and velocity vector v (evantually a
% different value of mu) and returns keplerian elements
%
% Inputs:
% r     : position vector           (km)
% v     : velocity vector           (km)
% mu    : gravitational parameter   (km^3/s^2)
%
% Outputs:
% a         : semi-major axis                       (km)
% ecc       : eccentricity
% i         : inclination                           (deg)
% W         : right ascension of ascending node     (deg)
% w         : argument of perigee                   (deg)
% theta     : true anomaly                          (deg)
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version

if nargin == 2
    mu = 398600.433;
end
%input check:
[l1,l2] = size(v);
if (l1<l2)
    v = v';
end
[l1,l2] = size(r);
if (l1<l2)
    r = r';
end
%angular momuntum:
h = cross(r,v); 
%inclination:
i = acos(h(end)/norm(h));
%eccentricity vector:
e = (cross(v,h)./mu)-(r/norm(r));
%eccentricity:
ecc= norm(e);
%specific energy:
E = 0.5*(v'*v)-(mu/norm(r));
%semi_major axis:
a = -mu./(2*E);

k = [0 0 1]; %k versor;
%Node line:
N = cross(k,h);
%RAAN:
if N(2)>=0
W = acos(N(1)/norm(N));
else
    W = 2*pi - acos(N(1)/norm(N));
end
%Argument of perigee:
if e(3)>=0
    w = acos(dot(N,e)/(norm(N)*ecc));
else
     w = 2*pi - acos(dot(N,e)/(norm(N)*ecc));
end
%radial velocity:
vr = dot(r,v)/norm(r);
%true anomaly:
if vr>=0
    theta = acos(dot (e,r)/(norm(r)*ecc));
else
    theta = 2*pi - acos(dot (e,r)/(norm(r)*ecc));
end