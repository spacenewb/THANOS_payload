function[Lx, Lz, pho_h, Ptx, delta_psi, delta_theta] = Ant_design(pho_x, pho_gr, g_swath, H, theta0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHYSICAL CONSTANTS
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIVEN DATA
f0 = 13.78e9; % Carrier frequency [Hz]
%theta0 = deg2rad(20); % Elevation Angle [rad]
%pho_x = 2.5; % Azimuth Resolution [m]
%pho_gr = 5; % Ground range Resolution [m]
sigma0_NESZ = 10^(0.1*(-22)); % Noise Equivalent Sigma Zero [dB]
%g_range_swath_max = 40e3; % Maximum Ground Swath [m]
F = 10^(0.1*(5)); % Reciever Noise Figure [dB]
sigma_0 = 10^(0.1*(-15)); % Sigma Zero [dB]
eta = 0.7; % Antenna Efficiency
duty_cycle = 0.2; % Duty Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSUMPTIONS
%H = 900e3; % orbit height [m]
v_orb = (Mu_B/(R_B+H))^0.5; % orbit velocity [m/s]
%g_swath = 10e3; % ground swath [m] (40 km maximum as per requirement)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACQUISITION AND FOCUSING OF SAR DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RADAR PARAMETERS
lambda = c/f0; % wavelength [m]
pho_r = pho_gr*sin(theta0); % range resolution [m]
pho_h = pho_gr*cos(theta0); % Alt resolution [m]
Lx = 2*pho_x; % physical antenna length [m]
delta_psi = lambda/Lx; % physical azimuth beamwidth
delta_theta = 2*asin(g_swath/(2*H*cos(theta0))); % calculate from ground swath [rad]
Lz = lambda/delta_theta; % physical antenna width [m] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLIGHT GEOMETRY
% ground range min and max
y_med = H*tan(theta0);
y_min = y_med - g_swath/2; % 1000? correct it
y_max = y_med + g_swath/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHETIC APERTURE
r_min = sqrt(y_min^2 + H^2); % slant range min
r_max = sqrt(y_max^2 + H^2); % slant range max
r = sqrt(y_med^2 + H^2); % slant range med
As = lambda/Lx*r_max; % maximum synthetic aperture
dx_ant = pho_x*.5; % spatial sampling of the synthetic aperture (0 < factor < 1)
% FOR SIMULATION
Ltot = 3*As; % Total length of the simulated orbit
xa = (-Ltot/2:dx_ant:Ltot/2); % Radar position along the orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL
PRI = dx_ant/v_orb; % Pulse Repetition Interval [s]
N_tau = As/(v_orb*PRI); % Number of Pulses in the Synthetic Aperture
Tg = duty_cycle*PRI; % Duration of each Pulse [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOISE
Ta_ = 290; % Antenna Physical Temperature [K]
T0 = 290; % Ambient Temperature [K]
T0_ref = 290; % Reference Temperature [K]
T_rx = (F-1)*T0_ref; % Reciever Temperature [K]
T_ant = (eta*Ta_ + (1-eta)*T0); % Antenna Temperature [K]
Tsys = T_ant + T_rx; % System Temperature [K]
N0 = (K*Tsys); % Noise Spectral Density 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANTENNA
theta_focal = theta0; % Antenna Focus Theta
psi_focal = 0; % Antenna Focus Psi
G_ant = (eta*4*pi/(delta_theta*delta_psi)); % Antenna Gain
G_ant_db = 10*log10(G_ant); % Antenna Gain [dB]
A_e = (G_ant*(lambda^2))/4/pi; % Effective Antenna Area [m^2]
% TRANSMITER
f_tx = (sinc((theta_focal-theta0)/delta_theta)*sinc(psi_focal/delta_psi))^2; % Directivity Function - Transmitter
Ptx = (sin(theta0)*N0*(4*pi*r_max^2)^2)/(sigma0_NESZ*pho_x*pho_r*N_tau*Tg*G_ant*A_e*(f_tx^2)); % Transmission Power

end
