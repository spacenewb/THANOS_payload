function [par] = tracks_fun_Titan(r,tspan,par)
% Returns coordinates (alpha, delta, latitude, longitude) and thetaG
%
% Inputs:
% r         : matrix position vectors
% t_span    : vector with intervals of time of evaluation   (sec)
% date      : initial date of evaluation                    (days or date)
% par       : structure with data (for constants)
%
% Outputs:
% par : structure that contains all the elements evaluated (alpha, delta,
% latitude, longitude, thetaG)
%
% Functions required:
% - CarToSph()
% - LongLat_date_fun()
% 
% Contributors:
% Gaballo Paolo
% 
% Versions:
% 2021-05-02, second version (Titan)
% 2021-02-13, first version (Earth)

% Conversion of r to right ascension ALPHA and declination DELTA
[delta,alpha] = CarToSph(r);            %[rad]

% Conversion to longitude LAMBDA and latitude PHI
% thetaG0 = 0;
T0 = par.const.T0;
[phi,lambda,thetaG_pert] = LongLat_date_fun(delta,alpha,tspan,T0);   %[rad]

% longitude and latitude in degrees
latitude = rad2deg(phi);                %[deg]
long = rad2deg(lambda);                 %[deg]

% longitute measured [-180,180] 
longitude= zeros(length(long),1);

for i=1:length(long)
    if long(i)>180
        longitude(i) = long(i) - 360;
    else
        longitude(i) = long(i);
    end
end

% Save data in structure orbit.coordinates
par.coordinates.alpha = alpha;
par.coordinates.delta = delta;
par.coordinates.longitude = longitude;
par.coordinates.latitude = latitude;
par.coordinates.thetaG = thetaG_pert;
end

