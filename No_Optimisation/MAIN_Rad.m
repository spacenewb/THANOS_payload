clc
clear all
close all

%% PHYSICAL CONSTANTS
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant

%% Media
f0 = 13.78e9; % Carrier frequency [Hz]
% ant_density = 4.2; % kg/m^2
% ref_density = 1.7; % kg/m^2 % gauge 1.1mm

%% Resolutions
ant_length_lim = 5; % [m]
rx_tx_ant = 1; % Ratio of rx antenna and Tx antenna length
% rx antenna longer to boost gain and hence SNR of rx
pho_r = 45; % Range resolution [m]
pho_h = 43; % Altitude resolution [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_antenna_length_tx = ant_length_lim/rx_tx_ant; % length of Tx antenna [m]

pho_x = max_antenna_length_tx/2; % Azimuth resolution [m]
theta_0 = acos(pho_h/pho_r); % Elevation Look Angle [rad]

pho_gr = pho_r*sin(theta_0); % Ground Range resolution [m]

%% DeltaV
Mass_SC = 6e3; % launch, Cassini dry:2523 [kg]
Expected_dry_mass = 2e3; %Our Satellite [kg]
Huygens_mass = 350;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g0 = 9.81; % [m/s^2]
Isp = 300; % [1/s]
DV = Isp*g0*log(Mass_SC/Expected_dry_mass); % [m/s]

Atmos_Alt = 900e3; % [m]
H_min = max([DeltaV_Alt_Incl(DV), Atmos_Alt]); % [m]
orb_T = (2*pi*(((R_B+H_min)^3)/Mu_B)^0.5); % [s]
v_orb = (Mu_B/(R_B+H_min))^0.5; % orbit velocity [m/s]

%% Mission Time
yr2sec = 31557600; % Years to seconds conversion factor

Acq_Modes = 5;
Mission_LifeTime_yrs = 8;

SAR_mission_time = 4.5*yr2sec;
Scatter_mission_time = 2*yr2sec;
Radio_mission_time = Scatter_mission_time*0;
CalibTest_mission_time = 0.3*yr2sec;
%Misc_mission_time = 1*yr2sec;

mission_time = SAR_mission_time+Scatter_mission_time+Radio_mission_time+CalibTest_mission_time; %+Misc_mission_time; % [s]
max_mission_extension = ((Mission_LifeTime_yrs*yr2sec)-mission_time); % End of life cap: 8 yrs

Primary_Mission_time_yrs = (SAR_mission_time+CalibTest_mission_time)/yr2sec;
Secondary_Mission_time_yrs = (Scatter_mission_time+Radio_mission_time)/yr2sec;

Mission_time_yrs = mission_time/yr2sec;
Max_mission_extension_yrs = max_mission_extension/yr2sec;

SAR_phases = 3;
Scatter_phases = 4; % One Before and after each Orb Inc Change
Radio_phases = 4; % One Before and after each Orb Inc Change
CalibTest_phases = 7; % One Before and after each Orb Inc Change, one during change of look angle

num_phases = SAR_phases+Scatter_phases+Radio_phases+CalibTest_phases;

%% Calculations for RADIOMETRY Mode only
Tx_DC = 0; % Duty Cycle of Tx Antenna
instr_DC = 1; % Payload Scan Orbital Duty Cycle
tx_orb_nos = 6; % number of orbits to transmit data of 1 orbit
transmission_window = 0.25; % Percentage of orbit time available for Tx
tx_SNR_min = 7; % Min. acceptable SNR for Tx [dB]
presum_factor = 1;
Processed_bandwidth = 1;

%% Compression and Encoding
%num_looks = 2;

% omegaT = 2*pi/(15.945421*24*3600); %[rad/s]
% dlambda = orb_T*omegaT;

num_phase_orbit = (Scatter_mission_time/Radio_phases/orb_T)
swath_angle = ((tx_orb_nos / instr_DC)*(2*pi))/(num_phase_orbit)
g_swath = round(tan(swath_angle/2)*2*R_B,3)
%g_swath = 2*pi*R_B/n_orbit_fullcoverage; % Ground Swath [m]

% plot
timing_acq = zeros(tx_orb_nos,1);
timing_acq(1) = instr_DC;

timing_obdp = zeros(tx_orb_nos,1);
timing_obdp(1) = (1-instr_DC);

timing_transmit  = ones(tx_orb_nos,1).*transmission_window;
timing_transmit(1) = 0;

timing = horzcat(horzcat(timing_acq, timing_obdp), timing_transmit);
% 
fig11 = figure;
bar(timing,0.5, 'stacked')
title('Timing Chart')
xlabel('Orbit Number');
ylabel('Duty Cycle');
legend('Scanning', 'On-Board Processing', 'Transmitting')
saveas(fig11, 'Graphs\2D\Timing.fig')
close(fig11)

%% Surface Coverage
orb_incl = deg2rad(45);
hh = R_B*(1 - sin(orb_incl));
hh_area = 2*pi*R_B*hh;
body_area = (4*pi*(R_B^2));
surface_coverage_hires = body_area - (2*hh_area);
surface_coverage_hires_percent = surface_coverage_hires/body_area;
surface_coverage_lowres_percent = 100;

%% Want to run script fast? Comment the next line.......
%tx_orb_nos = Tx_Lookup(tx_orb_nos, transmission_window, H_min, tx_SNR_min, Tx_DC, instr_DC, presum_factor, g_swath, 1.9e-5);

%%%%%%%%%

%% Reflector Design
ref_f = 2.5; % Focal length [m]
f_ant_off= 0.65;
numel_array = 4;
ref_B = 10e6;
Dz_factor = 3; % Ratio of required reflector width and actual width
Dz_offset_strut = 0.1;
Dx_factor = 1;

eta_a = 0.8; % Antenna Efficiency

%% SAR Design
[Lx, Lz, Dx, Dz, ref_depth, pho_h, pho_x, Ptx, Ptx_peak, Prx, ~, ApertureBlock, delta_psi, delta_theta, G_ant_dB_tx, G_ant_dB_rx, G_ref_dB, SNR_rx_dB, SNR_rc_dB, alfa, ant_BW, PRI, eta, duty_cycle, F_dB, sigma0_NESZ_dB, As, P_Rad, P_Rad_dB, SNR_Rad_dB] = Ant_design(Processed_bandwidth, presum_factor, f0, eta_a, pho_x, pho_gr, g_swath, H_min, theta_0, rx_tx_ant, ref_f, f_ant_off, numel_array);

Dx_final = Dx*Dx_factor;
Dz_final = (Dz*Dz_factor) + Dz_offset_strut;
Ant_Area = Lx*Lz;
Ant_res = min([pho_r, pho_h, pho_gr]);
Ant_Mass = Ant_Mass_eval(Acq_Modes, Ant_Area, Ant_res, 1);
Ref_Mass = Ant_Mass_eval(Acq_Modes, (pi*Dx*Dz_final), Ant_res, 1);

%% Beam Patterns and Shit

% Direction
theta_pointing = rad2deg(theta_0); % for S/C [deg] Nadir=0

% Incident direction
teta_i = 0; % Range Incidence[deg]
psi_i = 0.25; % Azimuth Incidence[deg]
array_spacing = 0.5; % Array spacing in terms of lambda

[max_teta_steering,max_psi_steering,N_z,N_x] = Beam_Patterns(Ptx,Ptx_peak,As,H_min,v_orb,theta_pointing,rad2deg(theta_0),teta_i,psi_i,array_spacing,f0,Lz,Lx);

Prx = ((1-Tx_DC)/Tx_DC)*(Prx*(numel_array^2))/(N_z);
Prx_peak = Prx/(1-Tx_DC);

%% Saving Antenna Struct

Ant.Array.Spacing = array_spacing;
Ant.Array.Numel_Range = N_z;
Ant.Array.Numel_Azimuth = N_x;
writetable(struct2table(Ant.Array), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Array');

Ant.Antenna.Length = Lx;
Ant.Antenna.Width = Lz;
Ant.Antenna.RecievedRadiometryPower = P_Rad;
Ant.Antenna.Efficiency = eta;
Ant.Antenna.CarrierFrequency = f0;
Ant.Antenna.Tx_Gain = G_ant_dB_tx;
Ant.Antenna.Rx_Gain = G_ant_dB_rx;
Ant.Antenna.BandWidth = ant_BW;
writetable(struct2table(Ant.Antenna), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Antenna');

Ant.Reflector.ApertureLength = Dx;
Ant.Reflector.ApertureWidth = Dz;
Ant.Reflector.Length = Dx_final;
Ant.Reflector.Width = Dz_final;
Ant.Reflector.Depth = ref_depth;
Ant.Reflector.FocalLength = ref_f;
Ant.Reflector.FocalOffset = f_ant_off;
Ant.Reflector.AppertureEfficiency = eta_a;
Ant.Reflector.FrequencyBandwidth = ref_B;
Ant.Reflector.GainMax = G_ref_dB;
Ant.Reflector.ApertureBlockage = ApertureBlock;
writetable(struct2table(Ant.Reflector), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Reflector');

Ant.Physical.Ant_Mass = Ant_Mass;
Ant.Physical.Ref_Mass = Ref_Mass;
writetable(struct2table(Ant.Physical), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Physical');

Ant.Scan.GroundSwath = g_swath;
Ant.Scan.OrbitAltitude = H_min;
Ant.Scan.OrbitVelocity = v_orb;
Ant.Scan.LookAngle = theta_0;
Ant.Scan.SyntheticApperture = As;
Ant.Scan.DutyCycle = instr_DC;
Ant.Scan.TimeOnOrbit = instr_DC*orb_T;
Ant.Scan.Processed_bandwidth = Processed_bandwidth;
Ant.Scan.Presum_factor = presum_factor;
writetable(struct2table(Ant.Scan), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Scan');

Ant.Res.Altitude = pho_h;
Ant.Res.Azimuth = pho_x;
Ant.Res.GroundRange = pho_gr;
Ant.Res.Range = pho_r;
writetable(struct2table(Ant.Res), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Resolution');

Ant.BW.Elevation = rad2deg(delta_theta);
Ant.BW.Azimuth = rad2deg(delta_psi);
writetable(struct2table(Ant.BW), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'BandWidth');

Ant.Pulse.PRI = PRI;
Ant.Pulse.DutyCycle = duty_cycle;
Ant.Pulse.ChirpRate = alfa;
writetable(struct2table(Ant.Pulse), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Pulse');

Ant.Noise.NoiseFigure_dB = F_dB;
Ant.Noise.NESZ_dB = sigma0_NESZ_dB;
Ant.Noise.SNR_Rad_dB = SNR_Rad_dB;
Ant.Noise.SNR_raw_dB = SNR_rx_dB;
Ant.Noise.SNR_RangeCompressed_dB = SNR_rc_dB;
writetable(struct2table(Ant.Noise), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Noise');

Ant.Tx.MinTransmissionSNR_dB = tx_SNR_min;
Ant.Tx.TransmissionWindow = transmission_window;
Ant.Tx.TransmissionOrbits = tx_orb_nos;
writetable(struct2table(Ant.Tx), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'DataTransmission');

Ant.Mission.MissionTime = Mission_time_yrs;
Ant.Mission.MissionExtension = Max_mission_extension_yrs;
Ant.Mission.EndOfLife = Mission_LifeTime_yrs;
Ant.Mission.MissionPhases = num_phases;
Ant.Mission.SAR_Mission_Time = SAR_mission_time;
Ant.Mission.Scatterometry_mission_time = Scatter_mission_time;
Ant.Mission.Radiometry_mission_time = Radio_mission_time;
Ant.Mission.CalibTest_mission_time = CalibTest_mission_time;
Ant.Mission.Primary_Mission_time_yrs = Primary_Mission_time_yrs;
Ant.Mission.Secondary_Mission_time_yrs = Secondary_Mission_time_yrs;
Ant.Mission.SAR_phases = SAR_phases;
Ant.Mission.Scatterometry_phases = Scatter_phases;
Ant.Mission.Radiometry_phases = Radio_phases;
Ant.Mission.CalibTest_phases = CalibTest_phases;
Ant.Mission.SurfaceCoverageHiRes = surface_coverage_hires_percent;
Ant.Mission.SurfaceCoverageLowRes = surface_coverage_lowres_percent;
writetable(struct2table(Ant.Mission), 'Rad_Design.xlsx', 'WriteMode', 'overwritesheet', 'Sheet', 'Mission');