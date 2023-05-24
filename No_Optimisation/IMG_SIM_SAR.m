function [Signal2Ambiguity, Duration_acquisition_simulated] = IMG_SIM_SAR(f0, orb_H, Orb_V, PRF, Tx_DC, Ptx_peak, Gtx, Grx, sigma0, L, theta0, delta_theta, pho_r, num_Targets)

%%% Orbital Velocity Amplified to show range migration
Orb_vel_mag_factor = 100;
Orb_V = Orb_vel_mag_factor*Orb_V;

%f0 = 13.78e9;
%orb_H = 900e3;
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
K = physconst('boltzmann'); % Boltzmann constant
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
%Orb_V = Orb_vel_mag_factor*(Mu_B/(R_B+orb_H))^0.5; % orbit velocity [m/s]

num_pulses = 3000;
dt = 1e-6;
%PRF = 3e3;
%Tx_DC = 0.2;

Duration_acquisition_simulated = num_pulses/PRF;
Synthetic_Aperture_Acquired_simulated = Duration_acquisition_simulated*Orb_V/Orb_vel_mag_factor;

%%
% Ptx_peak = 100/Tx_DC;
%Gtx = 10^(45/10);
%Grx = Gtx;
%sigma0 = 10^(-3/10);
%L = 10^(2/10);
lambda = c/f0;
%theta0 = deg2rad(32);
%delta_theta = deg2rad(4);
max_surf_alt = 5e3;

PRI = 1/PRF;
Tx_I = PRI*Tx_DC;
Tx_timevec = -PRI:dt:PRI;
Tx_amp = sinc((Tx_timevec)/(Tx_I/2));
%Tx_amp = Tx_amp*hamming(numel(Tx_amp));
Tx_signal = timeseries(Tx_amp, Tx_timevec);

X_interval = [(-((num_pulses-1)/2)*PRI*Orb_V) (+((num_pulses-1)/2)*PRI*Orb_V)];
Y_interval = [0 max_surf_alt];
Z_interval = [tan(theta0-(delta_theta/2)) tan(theta0+(delta_theta/2))].*orb_H;

R_interval = ((Z_interval.^2) + orb_H^2).^0.5; % Range limits of scanning

SAR_pos0 = [(-((num_pulses-1)/2)*PRI*Orb_V) orb_H 0];

%num_Targets = 10;
X_tgts = (rand(1,num_Targets).*abs(max(X_interval) - min(X_interval))) + min(X_interval);
Y_tgts = (rand(1,num_Targets).*abs(max(Y_interval) - min(Y_interval))) + min(Y_interval);
Z_tgts = (rand(1,num_Targets).*abs(max(Z_interval) - min(Z_interval))) + min(Z_interval);
Target_pos0 = [X_tgts', Y_tgts', Z_tgts'];

num_r_divs = (abs(max(Z_interval) - min(Z_interval))/pho_r);

%sampling_interval = 0.1*dt;
%Target_Times = ((2*R_interval(1)/c)-PRI):sampling_interval:((2*R_interval(2)/c)+PRI);
%Target_Times = ((2*R_interval(1)/c)):sampling_interval:((2*R_interval(2)/c));
Target_Times = linspace((2*R_interval(1)/c)-(Tx_I/2), (2*R_interval(2)/c)+(Tx_I/2), num_r_divs);
delta = zeros(size(Target_Times));

Img = zeros(num_pulses, numel(Target_Times));
Imgfft = zeros(num_pulses, numel(Target_Times));
for ii=1:num_pulses
    SAR_pos = SAR_pos0 + [((ii-1)*PRI*Orb_V) 0 0];
    for i=1:num_Targets
        TDist = sum((SAR_pos-Target_pos0(i,:)).^2, 'all')^0.5;
        Prx = (Ptx_peak*Gtx*Grx*(lambda^2)*sigma0)/(((4*pi)^2)*(TDist^2)*L);
        TDelay = (2*TDist/c);
        [~, minidx] = min(abs(Target_Times - TDelay));
        delta(minidx) = Prx;
    end
    resp = conv(delta, Tx_amp, 'same');
    respfft = conv(delta, abs(fftshift(fft(Tx_amp))), 'same');
    delta(:) = 0;
    Img(ii,:) = resp;
    Imgfft(ii,:) = respfft;
end

Img_hammed = zeros(size(Img));
h_window = hamming(floor(Tx_I*2/dt));
for hh = 1:num_pulses
    Img_hammed(hh,:) = conv( Img(hh,:), h_window, 'same' );
end
Tx_amp_hammed = conv( Tx_amp, h_window, 'same' );
Hamming_fft = abs(fftshift(fft(hamming(round(Tx_I/dt, 0)), 1000)));
LPA_BW = max(Hamming_fft(Hamming_fft>1)) - min(Hamming_fft(Hamming_fft>1));


%%
figure(1)
scatter(Target_pos0(:,3)./1e3, Target_pos0(:,1)./1e3, [], Target_pos0(:,2)./1e3, 'filled')
title('Simulated Ground Targets')
set(gca, 'YDir','reverse')
colormap(gray);
xlim(Z_interval./1e3)
ylim(X_interval./1e3)
xlabel('Range [km]')
ylabel('Azimuth [km]')
cb1 = colorbar;
ylabel(cb1, 'Target Altitude [km]')


figure(2)
sgtitle('Signal visualisation at various steps before processing for transmission')

subplot(2,3,1)
imshow(mat2gray(Img))
% imagesc(Img)
% colormap(gray)
title('Raw Image of Targets')

subplot(2,3,2)
imshow(mat2gray(Img_hammed))
%imagesc(Img_hammed)
title('Filtered Image of Targets')
colormap(gray)

subplot(2,3,3)
imshow(mat2gray(Imgfft))
%imagesc(Imgfft)
title('Processed (Uncompressed & Unencoded) Image of Targets')
colormap(gray)

subplot(2,3,4)
plot(Tx_signal,'k','LineWidth', 2);
title('Raw Pulse Response of Target @ Rx')
dim = [0.25 0 0 0.425];
ann = annotation('textbox',dim,'String',strcat('Raw Pulse Response'),'FitBoxToText','on');
ann.Color = 'red';
ann.FontSize = 8;
SAR_raw_dB = 10*log10(abs(trapz(Tx_amp(Tx_amp>0))/trapz(Tx_amp(Tx_amp<=0))));

subplot(2,3,5)
Tx_amp_hammed_norm = Tx_amp_hammed./max(Tx_amp_hammed);
Tx_signal_hammed_norm = timeseries(Tx_amp_hammed_norm, Tx_timevec);
plot(Tx_signal_hammed_norm,'k','LineWidth', 2);
title('Filtered Pulse Response of Target @ Rx')
dim = [0.525 0 0 0.425];
ann = annotation('textbox',dim,'String',strcat('Band-Pass Filter: ', num2str(round(LPA_BW,2)), ' Hz'),'FitBoxToText','on');
ann.Color = 'red';
ann.FontSize = 8;
SAR_hammed_dB = 10*log10(abs(trapz(Tx_amp_hammed(Tx_amp_hammed>0))/trapz(Tx_amp_hammed(Tx_amp_hammed<=0))));

subplot(2,3,6)
Tx_amp_fft = abs(fftshift(fft(Tx_amp)))./(max(abs(fftshift(fft(Tx_amp)))));
plot(timeseries(Tx_amp_fft, Tx_timevec),'k','LineWidth', 2);
title('Processed Pulse Response of Target @ OBDH')
dim = [0.804 0 0 0.425];
ann = annotation('textbox',dim,'String',strcat('Conv( Dirac\Delta, FFT(Pulse) )'),'FitBoxToText','on');
ann.Color = 'red';
ann.FontSize = 8;
SAR_conv_dB = 10*log10(abs(trapz(Tx_amp_fft(Tx_amp_fft>0))/trapz(Tx_amp_fft(Tx_amp_fft<=0))));

Signal2Ambiguity = [SAR_raw_dB SAR_hammed_dB SAR_conv_dB];

end
