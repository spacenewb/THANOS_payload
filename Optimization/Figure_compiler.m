clc
clear all

%% Put all the figures and data in a folder called "All_figs" inside current folder
%% Name of the figures must be: "1.fig", "2.fig"......
%% Name of the data files must be same as figure files: "1.mat", "2.mat"...
%% run the script
%% Find the optimum point in the final figure.
%% Get the coordinates of that point.
%% open the file named "Points_data.csv" in excel
%% Search the values of coordinates in the excel file.
%% Get the row number of the searched points.
%% type in matlab command window: "values(row_number, :)

%% the output is the required optimum design
%% the coordinates are the optimum design characteristics



Xvals = [];
Yvals = [];
Zvals = [];
Cvals = [];
values = [];

for i=1:19
    
    disp(i);
    
    fig_name = strcat('All_figs/',num2str(i),'.fig');
    points_name = strcat('All_figs/',num2str(i),'.mat');
    
    fig = openfig(fig_name);
    D = get (gca, 'Children' ); % get the handle of the line object
    
    points = load(points_name);
    
    points = points.OptimParetoFront;
    
    values = vertcat(values, points);
    
    XData = get (D, 'XData' ); % get the x data
    YData = get (D, 'YData' ); 
    ZData = get (D, 'ZData' ); 
    CData = get (D, 'CData' ); 
    
    Xvals = horzcat(Xvals, XData);
    Yvals = horzcat(Yvals, YData);
    Zvals = horzcat(Zvals, ZData);
    Cvals = vertcat(Cvals, CData);
    
    close all
end

Xvals = Xvals';
Yvals = Yvals';
Zvals = Zvals';

save('fig_data.mat','Xvals','Yvals','Zvals','Cvals')
save('points_data.mat','values')
writematrix(values,'Points_data.xlsx')


%%
figure(1)
pointsize=30;
scatter3(Xvals, Yvals, Zvals, pointsize, Cvals, 'filled');
title({'Filtered Optimum Objective Variables'; 'Corresponding to Optimum Pareto Front'})
xlim(gca,[0 10]) % Pixel Area
ylim(gca,[0 10]) % Antena Area
%set (gca, 'Color' , 'k' )
xlabel('Voxel Volume [m^3]')
ylabel('Antenna Area [m^2]')
zlabel('Mission Time [Years] ((GlobalScans HiRes = 2))')
cb = colorbar('eastoutside');
ylabel(cb, '-Rx SNR [-dB]' )

%% select the coordinates from the figure
x_value = 3.3809; % approximated at the fourth digit after the comma
y_value = 2.9070;
% z_value = ;
coords = horzcat(Xvals,Yvals,Zvals,Cvals);
x = find(round(coords(:,1),4) == x_value);
coords_x = coords(x,:);
y = find(round(coords_x(:,2),4) == y_value);
coords_y = coords_x(y,:); 
opt_row = x(y);

opt_coords = coords(opt_row,:) % optimal coordinates selected from the figure
opt_values = values(opt_row,:) % optimum design