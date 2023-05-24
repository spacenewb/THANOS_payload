% clc
% clear all

% f0 = 13*1e9; % Central frequency [Hz] Ku Band 12.5-18GHz
% B = 10e6;
% 
% D = 6; % Diameter [m]
% f = 1; % Focal length [m]
% eta_a = 0.7; % Apperture efficiency
% 
% [depth, subtend_angle, G_maxs] = heflector_sizing(f0, B, eta_a, D, f, 0)

function[depth, subtend_angle, G_maxs] = Reflector_sizing(f0, B, eta_a, D, f, plot_opt)

%% Initialise
depth = (D^2)/(16*f); % Reflector depth [m]
subtend_angle = atan(D/2/f); % Angle subtended by the Apperture

%% Directivity
%f0s = [12.5 15.25 13.78 18].*1e9;
f0s = [(f0-(B/2)) f0 (f0+(B/2))];

lambdas = physconst ( 'LightSpeed' )./f0s; % Wavelength [m]
G_maxs = 10.*log10( eta_a.* ((pi*D./lambdas).^2) ); % Gain [dB]

%% Directivity per angle
if plot_opt == 1
    G = zeros(1000,numel(lambdas));
    G_max = zeros(numel(lambdas),1);
    BW = zeros(numel(lambdas),1);

    hold on
    for i=1:numel(lambdas)
        [theta, G(:,i), G_max(i), BW(i)] = directivity(D, lambdas(i), eta_a);
        plot(theta, G(:,i))
    end
    hold off
end

end

