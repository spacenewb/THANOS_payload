function[theta, G, G_max,BW] = directivity(D, lambda, eta_a)

theta = linspace(-1.0,1.0,1000);

G_max = 10.*log10( eta_a.* ((pi*D./lambda).^2) ); % Gain [dB]

x = (pi*D/lambda).*sin(deg2rad(theta));

bessel = besselj(1,x);

G = G_max.*((2*bessel./x).^2);

k = 70;
BW = k*lambda/D; % Degrees (Half Power BW)

end
