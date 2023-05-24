function plot_orbit_ev(a0,e0,i0,RAAN0,w0,TA0,r_vect,TA_vect)
% Plots the evolution of the orbit over a time period starting from
% initial keplerian elements and position vector evolution
%
% Inputs :
% a0            : semi-major axis                       (km)
% e0            : eccentricity
% i0            : inclination                           (deg)
% RAAN0         : right ascension of ascending node     (deg)
% w0            : argument of perigee                   (deg)
% TA0           : true anomaly                          (deg)
% r_vect        : position vector [km]
% TA_vect       : velocity vector [deg]
% 
% Plots:
% Orbit evolution
% 
% Functions required:
% plot_orbit()
%
% Contributors:
% Duma Francesco
% Gaballo Paolo
% Merola Pierpaolo
% Parducci Alessandro
% 
% Versions:
% 2021-02-13, first version

% Plot nominal initial orbit
plot_orbit(a0,e0,i0,RAAN0,w0,TA0)
hold on

% Define r_vect components
x = r_vect(:,1);
y = r_vect(:,2);
z = r_vect(:,3);

% Define color over the periods using true anomaly variation
np_color = zeros(size(r_vect,1),1);
k = 1;
TA_vect = TA_vect * pi/180;     %[rad]
theta = unwrap(TA_vect);        %[rad]
theta = theta * 180/pi;         %[deg]
for i=1:size(r_vect,1)
    if theta(i) < k*360 + theta(1) 
        np_color(i) = k;
    else
        k = k+1;
        np_color(i) = k;
    end
end

% Plot 1 orbit over 50 periods
idx = find(mod(np_color-1,50) ~= 0);
x(idx)= NaN;
plot3(x,y,z);
xlabel('[km]','FontSize',9)
ylabel('[km]','FontSize',9)
zlabel('[km]','FontSize',9)
% title('Orbit Evolution')

% Set color properties
col = np_color;
patch(x,y,z,col,'FaceColor','none','EdgeColor','interp');
c = colorbar;
c.Limits = [0; col(end)];
c.Label;
c.Label.String = 'Periods';
% c.Label.Position = [0.47,3050,0];
% c.Label.Rotation = 0;
c.Label.FontSize = 10;

end
