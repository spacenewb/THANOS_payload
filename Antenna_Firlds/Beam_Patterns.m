% 
% %% Inputs
% 
% P_tx = 97.4549796;
% P_tx_peak = 487.274898; 
% As = 4150.85123; % Length of Synthetic Aperture [m]
% Orb_alt = 900e3; % Orbital Altitude [m]
% Orb_V = 1.6074e+03; % Orbital Velocity [m]
% 
% % DIRECTION OF RADIATION
% % Direction
% theta_pointing = 12; % for S/C [deg] Nadir=0
% theta_0 = 17; % SAR Look Angle [deg]
% 
% % Incident direction
% teta_i = 1; % Range Incidence[deg]
% psi_i = 0.25; % Azimuth Incidence[deg]
% array_spacing = 0.5; % Array spacing in terms of lambda
% 
% 
% 
% % ARRAY PARAMETERS
% f0 = 13.78e9; % carrier frequency [Hz]
% 
% La_z = 0.216; % Array Range aperture (theta) [m]
% La_x = 1.31; % Array Azimuth aperture (psi) [m]
% 
% 

function [] = Beam_Patterns(P_tx,P_tx_peak,As,Orb_alt,Orb_V,theta_pointing,theta_0,teta_i,psi_i,array_spacing,f0,La_z,La_x)

%% Computations

c = physconst('LightSpeed'); % propagation velocity  [m/s]
lambda = c/f0; % wavelength [m]

n_z = ((2*La_z/(lambda*array_spacing))-1)/3;
N_z = round(n_z,0); % Number of elements

n_x = ((2*La_x/(lambda*array_spacing))-1)/3;
N_x = round(n_x,0); % Number of elements

delta_z = La_z/N_z; % Range spacing between neighbouring elements
zn = ((-(N_z-1)/2:(N_z-1)/2)*delta_z); % Range position of each element
an_z = hamming(numel(zn)); % Coefficients of each Range elements
delta_x = La_x/N_x; % Azimuth spacing between neighbouring elements
xn = ((-(N_x-1)/2:(N_x-1)/2)*delta_x); % Azimuth position of each element
an_x = hamming(numel(xn)); % Coefficients of each Azimuth elements

% DIRECTION SCANNING
max_teta_steering = asind(lambda/2/La_z); % [deg]
max_psi_steering = asind(lambda/2/La_x); % [deg]

snz = zeros(numel(zn), numel(xn));
an = (an_x*an_z')';
for ii = 1:numel(zn)
    for jj = 1:numel(xn)
        sn_z = exp(-1i*2*pi/lambda*sind(teta_i)*zn(ii));
        sn_x = exp(-1i*2*pi/lambda*sind(psi_i)*xn(jj));
        % Combined signal: 
        snz(ii,jj) = an(ii,jj)*sn_z * sn_x;
    end
end


%%
% DIRECTION SCANNING
% Directions tested by setting the array coefficients
ang_res = 101; % How many points to evaluate per degree angle

teta_lim = round(5*max_teta_steering,0);
psi_lim = round(5*max_psi_steering,0);

num_angles_z = teta_lim*ang_res;
num_angles_x = psi_lim*ang_res;

teta = linspace(-teta_lim,teta_lim,num_angles_z);
psi = linspace(-psi_lim,psi_lim,num_angles_x);
% Corresponding spatial frequencies
fz = sind(teta)/lambda;
fx = sind(psi)/lambda;
% Coefficients
anz = exp(-1i*2*pi*zn(:)*fz);
anx = exp(-1i*2*pi*xn(:)*fx);
% Generation of the received signal
snzz = exp(-1i*2*pi/lambda*sind(teta_i)*zn(:));
snxx = exp(-1i*2*pi/lambda*sind(psi_i)*xn(:));
% Combined signal: 
% Fourier Transform of the signal received at each antenna 
sz = anz'*snzz;
sx = anx'*snxx;

s = sz*sx';

%% Plots
figure()
sgtitle(strcat(['Transmit Beam: [', num2str(teta_i), ', ', num2str(psi_i), '] Degrees'])) 
colormap(jet)

subplot(2,1,1)
imagesc(an)
colorbar
axis image
title('Array Elements Coefficients')
xlabel('Azimuth Elements')
ylabel('Range Elements')
xticks(linspace(1,N_x,5))
yticks(linspace(1,N_z,5))

subplot(2,1,2)
imagesc(atand(imag(snz)./abs(snz)))
colorbar
axis image
title('Signal Incidence Phase Angle - [deg]')
xlabel('Azimuth Elements')
ylabel('Range Elements')
xticks(linspace(1,N_x,5))
yticks(linspace(1,N_z,5))

figure()
sgtitle(strcat(['Array Response @ Focal Plane: [', num2str(teta_i), ', ', num2str(psi_i), '] Degrees'])) 
colormap(jet)

subplot(1,2,1)
imagesc(abs(s))
colorbar
axis image
title('Magnitude - Fourier Transpose')
xlabel('Azimuth Angle')
ylabel('Range Angle')
num_xticks = 5;
num_yticks = 5;
xticks((0:num_xticks).*((num_angles_x-1)/(num_xticks-1))+1)
xticklabels(linspace(psi(1),psi(end),num_xticks))
yticks((0:num_yticks).*((num_angles_z-1)/(num_yticks-1))+1)
yticklabels(linspace(teta(1),teta(end),num_yticks))

subplot(1,2,2)
imagesc(round(atand(imag(s)./abs(s)),2))
colorbar
axis image
title('Magnitude - Quadrature')
xlabel('Azimuth Angle')
ylabel('Range Angle')
xticks((0:num_xticks).*((num_angles_x-1)/(num_xticks-1))+1)
xticklabels(linspace(psi(1),psi(end),num_xticks))
yticks((0:num_yticks).*((num_angles_z-1)/(num_yticks-1))+1)
yticklabels(linspace(teta(1),teta(end),num_yticks))

figure()
plot(teta,abs(sz))
axis equal
title('Antenna Response - Range')
xlabel('Range Angle')
ylabel('Response Magnitude')

figure()
plot(psi,abs(sx))
axis equal
title('Antenna Response - Azimuth')
xlabel('Azimuth Angle')
ylabel('Response Magnitude')

%% Radiation Patterns
% DIRECTION SCANNING
% Directions tested by setting the array coefficients
ang_res = 5; % How many points to evaluate per degree angle

teta_lim = 180;
psi_lim = 180;

num_angles_z = teta_lim*ang_res;
num_angles_x = psi_lim*ang_res;

teta = linspace(-teta_lim,teta_lim,num_angles_z);
psi = linspace(-psi_lim,psi_lim,num_angles_x);
% Corresponding spatial frequencies
fz = sind(teta)/lambda;
fx = sind(psi)/lambda;
% Coefficients
anz = exp(-1i*2*pi*zn(:)*fz);
anx = exp(-1i*2*pi*xn(:)*fx);
% Generation of the received signal
snzz = exp(-1i*2*pi/lambda*sind(teta_i)*zn(:));
snxx = exp(-1i*2*pi/lambda*sind(psi_i)*xn(:));
% Combined signal: 
% Fourier Transform of the signal received at each antenna 
sz = anz'*snzz;
sx = anx'*snxx;

s = sz*sx';

%% Plots

figure()
sgtitle(strcat(['Transmit Beam: [', num2str(teta_i), ', ', num2str(psi_i), '] Degrees'])) 

subplot(2,1,1)
plot(teta,abs(sz))
title('Antenna Response - Range')
xlabel('Range Angle')
ylabel('Response Magnitude')

subplot(2,1,2)
plot(psi,abs(sx))
title('Antenna Response - Azimuth')
xlabel('Azimuth Angle')
ylabel('Response Magnitude')


figure()
polarplot(deg2rad(teta),abs(sz))
ax = gca;
ax.ThetaZeroLocation = 'top';
thetalim([-90 90])
title('Antenna Radiation Pattern - Range')


figure()
polarplot(deg2rad(psi),abs(sx))
ax = gca;
ax.ThetaZeroLocation = 'top';
thetalim([-90 90])
title('Antenna Radiation Pattern - Azimuth')


%% Transmission 

% Corresponding spatial frequency
fz = sind(theta_0+teta_i)/lambda;
% Coefficients
an = exp(-1i*2*pi*zn(:)*fz);
% FIELD CALCULATION
dy = 1e3; % spatial sampling along y [m]
dz = dy; % spatial sampling along z [m]
z_range = (-(Orb_alt/10):dz:Orb_alt); % vertical axis [m]
y_alt = (1:dy:Orb_alt); % horizontal axis [m]
[Z,Y] = ndgrid(z_range,y_alt); % Matrices of (z,y) coordinates

% field calculation
E_tx = P_tx_peak*As/Orb_V;
E = 0;
N = numel(zn);
for n = 1:N
    Rn = sqrt(Y.^2 + (Z-zn(n)).^2); % distances from the n-th antenna
    En = exp(-1i*2*pi/lambda*Rn)./Rn; % spherical wave from the n-th antenna
    E = E + an(n)*En; % total field
    E_factor = (cosd(theta_pointing)/cosd(theta_0));
    E = (E*E_factor);
end
E_db = 10*log10(abs(E.*E_tx));

%% Plots

figure(), 
imagesc(y_alt/1e3,z_range/1e3,E_db)
colormap(jet)
% shading interp
ax1 = gca;
ax1.YAxisLocation = 'right';
axis xy
hold on ;
line ([y_alt(1), 2*y_alt(end)], [0, 0], 'Color' , 'k', 'linewidth', 2);
line ([y_alt(1), 2*y_alt(end)*cosd(theta_pointing)], [0, 2*y_alt(end)*sind(theta_pointing)], 'Color' , 'r', 'linewidth', 3 );
line ([y_alt(1), 2*y_alt(end)*cosd(theta_0)], [0, 2*y_alt(end)*sind(theta_0)], 'Color' , 'k', 'linewidth', 1.5,'LineStyle','--' );
line ([y_alt(1), 2*y_alt(end)*cosd(theta_0+teta_i)], [0, 2*y_alt(end)*sind(theta_0+teta_i)], 'Color' , 'y', 'linewidth', 1,'LineStyle','--' );
hold off
camroll(-(90))
colorbar,
legend('NADIR Vector', 'S/C Pointing Vector', 'SAR Look Vector', 'Reciever Scanning Beam')
xlabel('Altitude [Km]'), 
ylabel('Range [Km]')
title('Radiated field intensity [dB]')

figure(), 
contourf(y_alt/1e3,z_range/1e3,E_db)
colormap(parula)
ax1 = gca;
ax1.YAxisLocation = 'right';
axis xy
camroll(-(90))
colorbar 
xlabel('Altitude [Km]'), 
ylabel('Range [Km]')
title('Radiated field pattern [dB]')

end
