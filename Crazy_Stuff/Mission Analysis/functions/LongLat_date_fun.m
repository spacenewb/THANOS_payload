function [phi,lambda,thetaG] = LongLat_date_fun(delta,alpha,tspan,T0)
% Computes latitude and longitude in rad from declination and right
% acension
% Considers initial date and intervals of time for the evaluation of the
% Greenwich sidereal time
%
% Inputs:
% delta     : declination                   (rad)
% alpha     : right ascension               (rad)
% date      : initial date, must be a vector [ y, mo, d, h, mi, s]
% tspan     : vector of intervals of time   (sec)
% 
% Outputs:
% phi       : latitude                      (rad)
% lambda    : longitude                     (rad)
% thetaG    : Greenwich sidereal time       (rad)
%
% Functions required:
% Greenwich_sid_t()
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version

%latitude
phi = delta;                                            
t0 = 0;                                             % initial time
omegaT = 2*pi/(15.945421*24*3600);   %[rad/s]                  % Earth's angular velocity

thetaG0 = 0;
thetaG = thetaG0 + omegaT .* (tspan-t0);            % Theta Greenwich
% [thetaG,~] = Greenwich_sid_t(date,tspan);             % Theta Greenwich
dim = size(thetaG);
dlambda = T0*omegaT
lambda=zeros(max(size(alpha)),1);

for i=1:length(tspan)
    
     lambda(i) = alpha(i) - thetaG(i);
     if lambda(i)<0
         lambda(i)=lambda(i)+2*pi*ceil(-lambda(i)/(2*pi));
     end

end
    
end
