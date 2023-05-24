function [tx_orb_nos] = Tx_Lookup(tx_orb_nos, transmission_window, Orb_alt, tx_SNR_min, Tx_DC, instr_DC,  presum_factor, g_swath, Tg)

BAQ_vec = [88 86 84 83 82];
%plot_var = 'datarate'; % 'time' or 'datarate'
%transmission_window = 0.5; % can be refined with LOS calculations
%tx_orb_nos = 2; % Number of orbits used for transmitting data of one orbit
% the "tx_orb_nos" value is also the mission-time multiplier factor
    
%Orb_alt = 1000e3;
%Tx_DC = 0.2;
%instr_DC = 0.2;

global data_rate_cap
global SNR_out_min
data_rate_cap = 4e6; % bits per second
SNR_out_min = tx_SNR_min; % dB

%% Initialising
pho = 5:10:200;
IN_snr = -10:1:39;

set(groot,'DefaultFigureColormap',jet)
prog_bar = waitbar((0/(length(IN_snr)*length(pho))),'Computing Values, Please wait...');
set(prog_bar,'WindowStyle','modal')
iter = 0;
iter_lim = (length(IN_snr)*length(pho)*numel(BAQ_vec));

Output_Comp_SNR = zeros(length(IN_snr), 2, length(BAQ_vec));
Output_Comp_SNR(1,:,:) = pho(1);
Output_SNR = zeros(length(IN_snr), length(pho), length(BAQ_vec));
Output_Tx = zeros(length(IN_snr), length(pho), length(BAQ_vec));
Output_SNR_frac = zeros(length(IN_snr), length(pho), length(BAQ_vec));
possibility = zeros(length(IN_snr), length(pho), length(BAQ_vec));

%% Computing
for i=1:numel(BAQ_vec)
    BAQ = BAQ_vec(i);
    
for j=1:length(IN_snr)
    IN_snr_dB = IN_snr(j);
    pho_r = pho(1);
    [~, b, ~] = data_performance(pho_r, presum_factor, g_swath, Tg, BAQ, IN_snr_dB, Orb_alt, Tx_DC, instr_DC, transmission_window, tx_orb_nos);
    Output_Comp_SNR(j,1,i) = b;
    pho_r = pho(end);
    [~, b, ~] = data_performance(pho_r, presum_factor, g_swath, Tg, BAQ, IN_snr_dB, Orb_alt, Tx_DC, instr_DC, transmission_window, tx_orb_nos);
    Output_Comp_SNR(j,2,i) = b; 
    
    for k=1:length(pho)
        pho_r = pho(k);
        [a, b, c] = data_performance(pho_r, presum_factor, g_swath, Tg, BAQ, IN_snr_dB, Orb_alt, Tx_DC, instr_DC, transmission_window, tx_orb_nos);
        Output_SNR(j,k,i) = b; 
        Output_Tx(j,k,i) = a/1e6; 
        possibility(j,k,i) = c;
        
        Output_SNR_frac(j,k,i) = abs(Output_SNR(j,k)/IN_snr_dB)*100;
        
        iter = iter+1;
        waitbar((iter/iter_lim));
    end
end




%% 2-D Plotting
xti = 2;
yti = 5;  

title_plt = strcat('Data Rate Req to Transmit Tx Win: ', num2str(transmission_window));
label_plt = '[Mbps]';

fig1 = figure('visible', 'off');

subplot(2,2,1)
pcolor(Output_SNR(:,:,i))
shading interp
title('BAQ SNR degradation')
xlabel('Range Resolution [m]');
ylabel('Antenna Reciever SNR [dB]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
cb = colorbar;
ylabel (cb, '[dB]' )

subplot(2,2,2)
pcolor(Output_Tx(:,:,i))
shading interp
title(title_plt)
xlabel('Range Resolution [m]');
ylabel('Antenna Reciever SNR [dB]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
cb = colorbar;
ylabel (cb, label_plt ) 

subplot(2,2,3)
pcolor(Output_SNR_frac(:,:,i))
shading interp
title('BAQ SNR degradation Percentage')
xlabel('Range Resolution [m]');
ylabel('Antenna Reciever SNR [dB]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
cb = colorbar;
ylabel (cb, '[% SNR_i_n]' )

subplot(2,2,4)
pcolor(Output_Comp_SNR(:,1:2,i))
shading interp
title4 = strcat('BAQ SNR degradation');
title(title4)
xlabel('BAQ Ratio');
ylabel('Antenna Reciever SNR [dB]');
BAQ_label_vec = ['8:8'; '8:6'; '8:4'; '8:3'; '8:2'];
xticks(1.5)
xticklabels(BAQ_label_vec(i,:))
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
cb = colorbar;
ylabel (cb, '[dB]' )

group_title = strcat('Data Performance @ Alt: ', num2str(Orb_alt/1e3), 'km');

sgtitle(group_title) 

fig_name = strcat('Graphs\2D\',num2str(BAQ),'.fig');

set(fig1, 'visible', 'on'); 
saveas(fig1, fig_name);
close(fig1)

end
close(prog_bar)
%% Lookup Table Calculation
max_BAQ_possible = zeros(length(IN_snr),length(pho));
a = zeros(length(BAQ_vec),1);
for i=1:length(IN_snr)
    for j=1:length(pho)
        for k=1:length(BAQ_vec)
            a(k) = possibility(i,j,k);
        end
        max_BAQ_possible(i,j) = max(a.*BAQ_vec');
    end
end

% Plotting
fig3 = figure;
pcolor(max_BAQ_possible)
%shading flat
colormap jet
caxis([80 90])
fig_title = strcat('FDBAQ Lookup Table [Tx Rate Lim:',num2str(data_rate_cap/1e6),'Mbps, Tx SNR Lim:',num2str(SNR_out_min),'dB]');
title(fig_title)
xlabel('Range Resolution [m]');
ylabel('Antenna Reciever SNR [dB]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
cb = colorbar;
cb.Ticks = [82, 83, 84, 86, 88];
cb.TickLabels = ['8:2'; '8:3'; '8:4'; '8:6'; '8:8']; 
ylabel (cb, 'BAQ Type' )
fig_name = 'Graphs\2D\FDBAQ_Lookup_Table';
saveas(fig3, fig_name);

%% Slice plot Trial

% fig2a = figure ;
% fs_a = slice(Output_SNR,[],[],1:numel(BAQ_vec)) ;
% shading interp
% set (gca, 'Color' , 'k' )
% set(fs_a,'FaceAlpha',0.8);
% title(strcat('BAQ SNR degradation @ Alt: ', num2str(Orb_alt/1e3), 'km'))
% xlabel('Range Resolution [m]');
% xticks(1:xti:numel(pho))
% xticklabels(pho(1:xti:numel(pho)))
% ylabel('Antenna Reciever SNR [dB]');
% yticks(1:yti:numel(IN_snr))
% yticklabels(IN_snr(1:yti:numel(IN_snr)))
% zlabel('BAQ Ratio');
% zticks(1:numel(BAQ_vec))
% zticklabels(BAQ_label_vec)
% cb = colorbar;
% ylabel (cb, '[dB]' )
% fig_name = strcat('Graphs\3D\BAQ_SNR_degradation.fig');
% saveas(fig2a, fig_name);

fig2b = figure ;
fs_b = slice(Output_SNR_frac,[],[],1:numel(BAQ_vec)) ;
shading interp
set (gca, 'Color' , 'w' )
set(fs_b,'FaceAlpha',0.8);
title(strcat('BAQ SNR degradation Percentage @ Alt: ', num2str(Orb_alt/1e3), 'km'))
xlabel('Range Resolution [m]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
ylabel('Antenna Reciever SNR [dB]');
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
zlabel('BAQ Ratio');
zticks(1:numel(BAQ_vec))
zticklabels(BAQ_label_vec)
cb = colorbar;
ylabel (cb, '[% SNR_i_n]' )
fig_name = strcat('Graphs\3D\BAQ_SNR_degradation_Percentage.fig');
saveas(fig2b, fig_name);

fig2c = figure ;
fs_c = slice(Output_Tx,[],[],1:numel(BAQ_vec)) ;
shading interp
set (gca, 'Color' , 'w' )
set(fs_c,'FaceAlpha',0.8);
title(strcat(title_plt, ' @ Alt: ', num2str(Orb_alt/1e3), 'km'))
xlabel('Range Resolution [m]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
ylabel('Antenna Reciever SNR [dB]');
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
zlabel('BAQ Ratio');
zticks(1:numel(BAQ_vec))
zticklabels(BAQ_label_vec)
cb = colorbar;
ylabel (cb, label_plt )
title_fig = 'Data_Transmit_Rate';
fig_name = strcat('Graphs\3D\', title_fig, '.fig');
saveas(fig2c, fig_name);

fig2d = figure ;
fs_d = slice(Output_Comp_SNR,[],[],1:numel(BAQ_vec)) ;
shading interp
set (gca, 'Color' , 'w' )
set(fs_d,'FaceAlpha',0.8);
title(strcat(title4, ' @ Alt: ', num2str(Orb_alt/1e3), 'km'))
xlabel('Range Resolution [m]');
xticks(1.5)
xticklabels('All Range Resolutions')
ylabel('Antenna Reciever SNR [dB]');
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
zlabel('BAQ Ratio');
zticks(1:numel(BAQ_vec))
zticklabels(BAQ_label_vec)
cb = colorbar;
ylabel (cb, '[dB]' )
fig_name = strcat('Graphs\3D\BAQ_SNR_Degradation.fig');
saveas(fig2d, fig_name);

fig2e = figure ;
fs_e = slice(possibility,[],[],1:numel(BAQ_vec)) ;
shading flat
set (gca, 'Color' , 'w' )
set(fs_e,'FaceAlpha',0.8);
title(strcat('Tx Possibility @ Alt: ', num2str(Orb_alt/1e3), 'km'))
xlabel('Range Resolution [m]');
xticks(1:xti:numel(pho))
xticklabels(pho(1:xti:numel(pho)))
ylabel('Antenna Reciever SNR [dB]');
yticks(1:yti:numel(IN_snr))
yticklabels(IN_snr(1:yti:numel(IN_snr)))
zlabel('BAQ Ratio');
zticks(1:numel(BAQ_vec))
zticklabels(BAQ_label_vec)
cb = colorbar;
ylabel (cb, '[Possibility of Data Transmission]' )
fig_name = strcat('Graphs\3D\Tx_Possibility.fig');
saveas(fig2e, fig_name);

set(groot,'DefaultFigureColormap',parula)

close all

end