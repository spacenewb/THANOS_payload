%% SCATTEROMETRY
clc
clear all
close all

% INPUT PARAMETERS
f0 = 1.378000000000000e+10  % Carrying Frequency [Hz]
c = physconst('lightspeed'); % wave propagation velocity
K = physconst('boltzmann'); % Boltzmann constant
lambda = c/f0
H = 900e3 % orbit height [m]
theta0 = deg2rad(0) % pointing angle in elevation [rad]
teta_point = theta0;
Le = 5;  % antenna height [m]
Lx = 3; % antenna width [m]
ant_eff = 0.8 % antenna efficiency
Ptx = 1; % transmitted power [W]

%Resolution (dx = dy)
d_x = lambda/Lx * H            %[m]
d_e = lambda/Le * H               %[m]

sigma0_NESZ_dB = -15;
sigma0_NESZ = 10^((sigma0_NESZ_dB)/10); % Noise Equivalent Sigma Zero [dB]
F_dB = 3 % Reciever Noise Figure [dB]
F = 10^(F_dB/10); % Reciever Noise Figure
%sigma0_dB = 10*log10(sigma0_eval(rad2deg(theta0))); %-3; % Nature Paper
sigma0_dB = -3 % Nature Paper
sigma_0 = 10^((sigma0_dB)/10); % Sigma Zero [dB]
%%
% sigma0_dB = -15; % backscatter coefficient
% backscatter value taken from:
% A. Freeman and S. L. Durden, "A three-component scattering model for
% polarimetric SAR data," in IEEE Transactions on Geoscience and Remote Sensing, 
% vol. 36, no. 3, pp. 963-973, May 1998, doi: 10.1109/36.673687
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUND FOOTPRINT
delta_teta = lambda/Le; % beamwidth in elevation [rad]
delta_psi = lambda/Lx; % beamwidth in azimuth [rad]
% max and min elevation angles (w.r.t. Nadir)
teta_max = teta_point + delta_teta/2;
teta_min = teta_point - delta_teta/2;
% Max and min distances
R_min = H/cos(teta_min);
R_max = H/cos(teta_max);
% Ground footprint extension
x_min = -R_max*delta_psi/2; 
x_max = +R_max*delta_psi/2; 
y_min = H*tan(teta_min); 
y_max = H*tan(teta_max); 
% Spatial axes (with some margin in both direction)
margin_x = 3e3;
x_ax = linspace(x_min-margin_x,x_max+margin_x,501);
margin_y = 1e4;
y_ax = linspace(y_min-margin_y,y_max+margin_y,501);
% Spatial sampling
dx = x_ax(2) - x_ax(1);
dy = y_ax(2) - y_ax(1);
% Spatial coordinates of ground points
[Y,X] = ndgrid(y_ax,x_ax);
% Range at each ground position
R = sqrt(X.^2 + Y.^2 + H.^2);
% Azimuth angle at each ground position
PSI = asin(X./R); % X = R*sin(psi)
% Elevation angle at each ground position
TETA = asin(Y./R./cos(PSI)); % Y = R*cos(psi)*sin(teta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANTENNA PATTERN - approximated as the product of two squared sinc
% Tx pattern
f_tx = (sinc((TETA-teta_point)/delta_teta).*sinc(PSI/delta_psi)).^2;
% Rx pattern - the same antenna is used in Tx and Rx
f_rx = f_tx; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RADAR EQUATION
A = Le*Lx; % antenna physical area
G = ant_eff*4*pi*A/(lambda^2); % antenna gain
Ae = (lambda^2)/4/pi*G; % antenna effective area
G_dB = 10*log10(G); % Gain in dB

B_x = c/2/d_x; % bandwidth [Hz]
B_e = c/2/d_e; % bandwidth [Hz]

B = max([B_x B_e])

% Incident power density at the surface
Si_scatt = Ptx./(4*pi*R.^2)*G.*f_tx;
Si = Si_scatt;
% Scattered power density at the receiver from each ground position with area dA = dx*dy
Ss = Si./(4*pi*R.^2)*sigma_0*dx*dy;
% Power at the receiver form each ground position with area dA = dx*dy
Prx_xy = Ss*Ae.*f_rx;
% Received power in Watts
Prx_scatt = sum(sum(Prx_xy));


 %% NOISE

Ta_ = 273; % Antenna Physical Temperature [K]
T0 = 273; % Ambient Temperature [K]
T0_ref = 95; % Reference Temperature [K] % Titan
T_rx = (F-1)*T0_ref; % Reciever Temperature [K]
T_ant = (ant_eff*Ta_ + (1-ant_eff)*T0); % Antenna Temperature [K]
Tsys = T_ant + T_rx; % System Temperature [K]
N0 = (K*Tsys); % Noise Spectral Density 
Pnoise = N0*B; % Noise Power 
P_noise_dB = 10*log10(Pnoise);

% Power at receiving Radiometry antenna
Prx_rad = K*T0_ref*B;          %[W]


% SNR
SNR_scatt = (Prx_scatt/Pnoise);
SNR_rad = (Prx_rad/Pnoise);

SNR_scatt_dB = 10*log10(SNR_scatt)
SNR_rad_dB = 10*log10(SNR_rad)

figure, imagesc(x_ax*1e-3,y_ax*1e-3,Si), axis xy 
colorbar
xlabel('x [Km]'), ylabel('y [Km]')
title('Incident power density at the surface [W/m^2]')

 

figure, imagesc(x_ax*1e-3,y_ax*1e-3,Prx_xy), axis xy 
colorbar
xlabel('x [Km]'), ylabel('y [Km]')
title('Received power from individual ground positions [W]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%