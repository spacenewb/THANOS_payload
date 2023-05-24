function[Lx_tot, Lz_tot, D_x, D_z, depth_max, pho_h, pho_x, Ptx, Ptx_peak, Prx, Prx_peak, Aperture_Block, delta_psi, delta_theta, G_db_tx, G_db_rx, G_ref_dB, SNR_rx_dB, SNR_rc_dB, alfa, B, PRI, eta_a, duty_cycle, F_dB, sigma0_NESZ_dB, As, P_Rad, P_Rad_dB, SNR_Rad_dB, PRF, sigma_0, L] = Ant_design(Processed_bandwidth, presum_factor, f0, eta_a, pho_x, pho_gr, g_swath, orb_H, theta0, rx_tx_ant, focal_l, f_ant_off, numel_array)


%% PHYSICAL CONSTANTS
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant

%% ASSUMPTIONS
%f0 = 13.78e9; % Carrier frequency [Hz]
sigma0_NESZ_dB = -15;
sigma0_NESZ = 10^((sigma0_NESZ_dB)/10); % Noise Equivalent Sigma Zero [dB]
F_dB = 3; % Reciever Noise Figure [dB]
F = 10^(F_dB/10); % Reciever Noise Figure
sigma_0_dB = 10*log10(sigma0_eval(rad2deg(theta0))); %-3; % Nature Paper
sigma_0 = 10^((sigma_0_dB)/10); % Sigma Zero [dB]
duty_cycle = 0.2; % Duty Cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_orb = (Mu_B/(R_B+orb_H))^0.5; % orbit velocity [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = 2*pho_x; % physical antenna length [m]
B_doppler = 2*v_orb/Lx;
%presum_factor = 0.4;
PRF = presum_factor*B_doppler;
degrad_factor = presum_factor*Processed_bandwidth;
Lx_eff = Lx/degrad_factor;
swathDegradFactor = Lx/Lx_eff;
pho_x = Lx_eff/2;
%SNR_S_A_dB = 5+(20*(presum_factor + 0.15)*Processed_bandwidth);
g_swath = g_swath*swathDegradFactor;

%% RADAR PARAMETERS
theta0 = abs(theta0);
lambda = c/f0; % wavelength [m]
pho_r = pho_gr/sin(theta0); % range resolution [m]
pho_h = pho_r*cos(theta0); % Alt resolution [m]
%Lx = 2*pho_x; % physical antenna length [m]
delta_psi = lambda/Lx; % physical azimuth beamwidth
delta_psi_rx = lambda/(Lx*rx_tx_ant); % physical azimuth beamwidth Rx

delta_theta_eq = @(del_theta) (g_swath - 2*orb_H*sin(del_theta)) / (cos(2*theta0) + cos(del_theta));
options_f = optimoptions ('fsolve', 'Display', 'none');
delta_theta = fsolve(delta_theta_eq,0,options_f);

Lz = lambda/delta_theta; % physical antenna width [m] 
L_dB_km = 0.015;
PRF_ul = c/2/g_swath/sin(theta0);
%% FLIGHT GEOMETRY
y_max = orb_H*tan(theta0+(delta_theta/2));
y_min = y_max - g_swath;

%% SYNTHETIC APERTURE
r_max = sqrt(y_max^2 + orb_H^2); % slant range max
r_min = sqrt(y_min^2 + orb_H^2); % slant range min
As = lambda/Lx*r_max; % maximum synthetic aperture
dx_ant = pho_x*0.5; % spatial sampling of the synthetic aperture (0 < factor < 1)
x_amb = lambda/2/dx_ant*r_min; % ambiguous along-track distance if no antenna pattern is used
% FOR SIMULATION
Ltot = 3*As; % Total length of the simulated orbit
xa = (-Ltot/2:dx_ant:Ltot/2); % Radar position along the orbit
Nx = length(xa);

%% SIGNAL

PRI = dx_ant/v_orb; % Pulse Repetition Interval [s]
PRI = 1/PRF;

PRI_min = 2*(r_max-r_min)/c; % To avoid Range Ambiguities [s]
PRI_max = Lx/2/v_orb; % To avoid Azimuth Ambiguities [s]
%PRI = (PRI_min + PRI_max)/2;

N_tau = As/(v_orb*PRI); % Number of Pulses in the Synthetic Aperture
Tg = duty_cycle*PRI; % Duration of each Pulse [s]

B = c/2/pho_r; % bandwidth [Hz]
alfa = B/Tg; % chirp rate [Hz^2]
SNR_dB = 10*log10(sigma_0/sigma0_NESZ); % Signal to noise ratio [dB]
SNR = sigma_0/sigma0_NESZ;

%% NOISE
Ta_ = 273; % Antenna Physical Temperature [K]
T0 = 273; % Ambient Temperature [K]
T0_ref = 95; % Reference Temperature [K] % Titan
T_rx = (F-1)*T0_ref; % Reciever Temperature [K]
T_ant = (eta_a*Ta_ + (1-eta_a)*T0); % Antenna Temperature [K]
Tsys = T_ant + T_rx; % System Temperature [K]
N0 = (K*Tsys); % Noise Spectral Density 
Pnoise = N0*B; % Noise Power 
P_noise_dB = 10*log10(Pnoise);

P_Rad = K*T0_ref*B
P_Rad_dB = 10*log10(Pnoise);

SNR_rad = P_Rad/Pnoise;
SNR_Rad_dB = 10*log10(SNR_rad)

%% ANTENNA AND REFLECTOR


%% Ant Ref Design


D_x = lambda/delta_psi;
D_z = lambda/delta_theta;

G_ref_dB = 10.*log10(((pi/delta_theta)*(pi/delta_psi))); % Gain [dB]

depth_max = (max([D_x, D_z])^2)/(16*focal_l); % Reflector depth [m]

delta_psi_ant = 2*atan(D_x/2/focal_l);
delta_theta_ant = 2*atan(D_z/2/focal_l);
G_ant_dB = 10.*log10(((pi/delta_theta_ant)*(pi/delta_psi_ant)) ); % Gain [dB]


Lx = lambda/delta_psi_ant;
Lz = lambda/delta_theta_ant;

scale_ang_x = atan(((D_x - Lx)/2)/focal_l);
scale_ang_z = atan(((D_z - Lz)/2)/focal_l);

array_delta_psi = 2*scale_ang_x;
array_delta_theta = 2*scale_ang_z;
G_array_dB = 10.*log10(((pi/array_delta_theta)*(pi/array_delta_psi)) ); % Gain [dB]
G_array = 10^(G_array_dB/10);

delta_Lx = tan(scale_ang_x)*f_ant_off*2;
delta_Lz = tan(scale_ang_z)*f_ant_off*2;

Lx_tot = Lx + delta_Lx;
Lz_tot = Lz + delta_Lz;

Ae_array = Lx_tot*Lz_tot;

block_percent = (Ae_array/(D_x*D_z))*100;

%%



L_dB = L_dB_km*(r_max/1e3);
L = 10^(L_dB/10);
RCS = sigma_0*r_max*pho_x*pho_r/sin(theta0);

% theta_focal = theta0; % Antenna Focus Theta
% psi_focal = 0; % Antenna Focus Psi
G_array_tx = G_array;
G_array_rx = G_array*(numel_array^2);

Ptx = (SNR*L*N0*(4*pi*r_max^2)^2) / ((RCS*Tg*G_array_tx*Ae_array)); % Transmission Power
Prx = (SNR*L*N0*(4*pi*r_max^2)^2) / ((RCS*Tg*G_array_rx*Ae_array)); % Transmission Power

Ptx_peak = Ptx/duty_cycle;
Prx_peak = Prx/duty_cycle;

%% SNR Predicted
Prxed = (Ptx*G_array_rx*Ae_array*RCS)/(L*(4*pi*r_max^2)^2);
% Before range compression
SNR_rx_dB = 10*log10(Prxed/Pnoise); % Terrible values due to very low Prx and slightly high Pnoise
% After range compression
SNR_rc_dB = 10*log10(Prxed*Tg/N0); % Terrible values due to very low Prx and slightly high Pnoise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fake outputs
Aperture_Block=block_percent/100;
gain_antenna_tx=G_array_tx;
gain_antenna_rx=G_array_rx;
G_db_tx=10*log10(gain_antenna_tx);
G_db_rx=10*log10(gain_antenna_rx);
end