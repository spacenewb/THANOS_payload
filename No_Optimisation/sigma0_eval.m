function[sigma_0] = sigma0_eval(querry_theta0)

data = load('sigma0theta_ref.mat');

theta_0s = table2array(data.Sigma0(:,1));
sigma_0s = table2array(data.Sigma0(:,2));


sigma_0 = interp1(theta_0s, sigma_0s, querry_theta0);

%% Plot

figc = figure;
plot(theta_0s, sigma_0s, 'LineWidth', 2, 'Color', 'r')
shading interp
title('Radar Backscattering Coefficient')
xlabel('Look Angle \theta_0 [Degrees]');
ylabel('Radar Backscattering Coefficient \sigma_0');
set (gca, 'YScale' , 'log' )
% xticks(1:xti:numel(pho))
% xticklabels(pho(1:xti:numel(pho)))
% yticks(1:yti:numel(IN_snr))
% yticklabels(IN_snr(1:yti:numel(IN_snr)))
% cb = colorbar;
% ylabel (cb, '[% SNR_i_n]' )
grid on

fig_name = strcat('Graphs\2D\Sigma0_Theta0.fig');

saveas(figc, fig_name);
close(figc)

end
