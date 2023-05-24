function [req_data_rate, SNR_degradation, possible] = data_performance(pho_r, presum_factor, g_swath, Tg, BAQ, IN_snr_dB, Orb_alt, Tx_DC, instr_DC, transmission_window, tx_orb_nos)

%% Declared Variables
Rx_channel = 1;
dir_channel = 2; % one for phi and one for theta
Tx_channel = 1;

Array_Rg = 3;
Array_Az = 3;

%sampling_f = 533.33e6;

quantisation = 8; % bits
direction_quantisation = 8; % bits

%% Physical Constants
%peak_data_rate = 430e3; % bps

body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant

v_orb = (Mu_B/(R_B+Orb_alt))^0.5; % orbit velocity [m/s]
Orb_period = (2*pi*(((R_B+Orb_alt)^3)/Mu_B)^0.5); % s

%% Raw Data Rate
Bandwidth = c/2/pho_r;
Nyquist_sampling_freq = 2*Bandwidth;
sampling_f = Nyquist_sampling_freq; % Hz

Array_tot = Array_Rg*Array_Az;
Rx_DC = 1 - Tx_DC;

raw_data_1channel = Rx_channel*Rx_DC*(sampling_f*quantisation)*presum_factor*((g_swath/c) + Tg); % bps
raw_data_allChannel = Array_tot*raw_data_1channel*presum_factor; % bps

raw_dir_data = dir_channel*Rx_DC*(sampling_f*direction_quantisation)*presum_factor; % bps

processed_data = raw_data_1channel + raw_dir_data;

%% Compression and Encoding - BAQ [8:8, 8:6, 8:4, 8:3, 8:2] Lloyd Max Algorithm
% Table-A: [SNR]; (x:SNR_in_dB, y:SNR_BAQ_dB) @ Gamma_clip=-8.8dB
% Xrange(-10:40)

% Table-B: [Clip]; (x:Gamma_clip_dB, y:BAQ_conversion_gain_dB) @ SNR_in=0dB
% Xrange(-14:-2) but for 82(-14:-4.8)

% input_type: SNR or Clip (to decide lookup table)
input_type = 'SNR';

[outputY, output_type, BAQ_Compression_ratio] = BAQ_LookupTable(BAQ, input_type, IN_snr_dB);

if output_type == 'BAQ_SNR'
        BAQ_snr_dB = outputY;
        BAQ_gain_dB = 'not evaluated';
        BAQ;
        SNR_degradation = (BAQ_snr_dB - IN_snr_dB);
elseif output_type == 'BAQ_Gain'
        BAQ_snr_dB = 'not evaluated';
        BAQ_gain_dB = outputY;
        BAQ_gain_dB;
end

%% Transmit (Orbit level) calculations
DATA_tx = BAQ_Compression_ratio*processed_data;

orb_data = DATA_tx*Orb_period*instr_DC; % bits of data generated in one orbit

%transmission_window = 0.5;
req_data_rate = orb_data/(transmission_window*tx_orb_nos*Orb_period);

global data_rate_cap
global SNR_out_min

if req_data_rate <= data_rate_cap && BAQ_snr_dB>=SNR_out_min
    possible = 1;
else
    possible = 0;
end

end