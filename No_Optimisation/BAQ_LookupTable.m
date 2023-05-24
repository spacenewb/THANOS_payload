function [outputY, output_type, Compression_ratio] = BAQ_LookupTable(BAQ, input_type, inputX)

% Assumptions: Signal is random variable normally distributed with p(x) = 1
% Noise is modelled as gaussian white noise
% Quantisation noise is considering only level error, not granulation

%% Decalare read data type

% Table A (x:SNR_in_dB, y:SNR_BAQ_dB) @ Gamma_clip=-8.8dB
% Table B (x:Gamma_clip_dB, y:BAQ_conversion_gain_dB) @ SNR_in=0dB

BAQ_type = num2str(BAQ);

switch input_type
    case 'SNR'
        table_type = 'A';
        output_type = 'BAQ_SNR';
    case 'Clip'
        table_type = 'B';
        output_type = 'BAQ_Gain';
end

table_directory = strcat('CSV\',BAQ_type,table_type,'.csv');

%% Read Data

Table = readmatrix(table_directory);
DataX = Table(:,1);
DataY = Table(:,2);

outputY = interp1(DataX, DataY, inputX);

%% Maintain BAQ_snr <= In_snr

if input_type == 'SNR'
    outputY = min(outputY, inputX);
end

%% Compression Ratio

Compression_ratio = str2num(BAQ_type(2))/str2num(BAQ_type(1));

end

