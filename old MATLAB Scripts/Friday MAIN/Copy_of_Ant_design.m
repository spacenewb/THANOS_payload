function[Lx, Lz, pho_h, G_ant_db, Ptx] = Ant_design(pho_x, pho_gr, g_swath)

%% Global Variable Calling
global c K
global f0 theta0 sigma0_NESZ F sigma_0 eta duty_cycle
global H v_orb

%% ACQUISITION AND FOCUSING OF SAR DATA

%% RADAR PARAMETERS
lambda = c/f0; % wavelength [m]
pho_r = pho_gr*sin(theta0); % range resolution [m]
pho_h = pho_gr*cos(theta0); % Alt resolution [m]
B = c/2/pho_r; % bandwidth [Hz]
Lx = 2*pho_x; % physical antenna length [m]
delta_psi = lambda/Lx; % physical azimuth beamwidth
delta_theta = 2*asin(g_swath/(2*H*cos(theta0))); % calculate from ground swath [rad]
Lz = lambda/delta_theta; % physical antenna width [m] 
AR = Lx/Lz;
voxel_size = [pho_x pho_gr pho_h];

%% FLIGHT GEOMETRY
% ground range min and max
y_med = H*tan(theta0);
y_min = y_med - g_swath/2; % 1000? correct it
y_max = y_med + g_swath/2;

%% SYNTHETIC APERTURE
r_min = sqrt(y_min^2 + H^2); % slant range min
r_max = sqrt(y_max^2 + H^2); % slant range max
r = sqrt(y_med^2 + H^2); % slant range med
As = lambda/Lx*r_max; % maximum synthetic aperture
dx_ant = pho_x*.5; % spatial sampling of the synthetic aperture (0 < factor < 1)
x_amb = lambda/2/dx_ant*r_min; % ambiguous along-track distance if no antenna pattern is used

%% FOR SIMULATION
Ltot = 3*As; % Total length of the simulated orbit
xa = (-Ltot/2:dx_ant:Ltot/2); % Radar position along the orbit
Nx = length(xa);

%% SIGNAL
PRI = dx_ant/v_orb; % Pulse Repetition Interval [s]
PRF = 1/PRI; % Pulse Repetition Frequency [Hz]

PRI_min = 2*(r_max-r_min)/c; % To avoid Range Ambiguities [s]
PRI_max = Lx/2/v_orb; % To avoid Azimuth Ambiguities [s]
PRF_max = 1/PRI_min; % To avoid Azimuth Ambiguities [Hz]
PRF_min = 1/PRI_max; % To avoid Range Ambiguities [Hz]

N_tau = As/(v_orb*PRI); % Number of Pulses in the Synthetic Aperture
Tg = duty_cycle*PRI; % Duration of each Pulse [s]
alfa = B/Tg; % chirp rate [Hz^2]
SNR = 10*log10(sigma_0/sigma0_NESZ); % Signal to noise ratio [dB]

%% NOISE
Ta_ = 290; % Antenna Physical Temperature [K]
T0 = 290; % Ambient Temperature [K]
T0_ref = 290; % Reference Temperature [K]
T_rx = (F-1)*T0_ref; % Reciever Temperature [K]
T_ant = (eta*Ta_ + (1-eta)*T0); % Antenna Temperature [K]
Tsys = T_ant + T_rx; % System Temperature [K]
N0 = (K*Tsys); % Noise Spectral Density 
Pnoise = N0*B; % Noise Power 

%% ANTENNA
theta_focal = theta0; % Antenna Focus Theta
psi_focal = 0; % Antenna Focus Psi
G_ant = (eta*4*pi/(delta_theta*delta_psi)); % Antenna Gain
G_ant_db = 10*log10(G_ant); % Antenna Gain [dB]
A_e = (G_ant*(lambda^2))/4/pi; % Effective Antenna Area [m^2]
% TRANSMITER
f_tx = (sinc((theta_focal-theta0)/delta_theta)*sinc(psi_focal/delta_psi))^2; % Directivity Function - Transmitter
Ptx = (sin(theta0)*N0*(4*pi*r_max^2)^2)/(sigma0_NESZ*pho_x*pho_r*N_tau*Tg*G_ant*A_e*(f_tx^2)); % Transmission Power
% RECIEVER
f_rx = f_tx; % Directivity Function - Reciever (Same Antenna)
Prx = Ptx*G_ant*A_e*(f_rx^2)*sigma_0*pho_x*pho_r/(sin(theta0)*(4*pi*r^2)^2); % Reciever Power

%% PREDICTED SNR
% Before range compression
SNR_rx_dB = 10*log10(Prx/Pnoise); % Terrible values due to very low Prx and slightly high Pnoise
% After range compression
SNR_rc_dB = 10*log10(Prx*Tg/N0); % Terrible values due to very low Prx and slightly high Pnoise

end
