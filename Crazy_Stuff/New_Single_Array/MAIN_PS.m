clc
clear all

%% Nomenclature
%DesignVariables: The variables that we use to design the SAR eg:Alt,swath,elevation..
%ObjectiveVariable: The variable used to judge the optimisation score eg:power,mass..
%ParetoFront: The set of Design variables tht are a result of Pareto optimum ObjectiveVariables

%% Create Parallel Pool
if gcp('nocreate')==0
parpool()
end

%% Limits
voxel_limits = [25 26];
Ant_a_limits = [11 12];

%% Fantasy DeltaV
Mass_Cassini = 6e3; % launch, dry:2523
Expected_dry_mass = 2e3; %Our Satellite

g0 = 9.81;
Isp = 300;
DV = Isp*g0*log(Mass_Cassini/Expected_dry_mass);

%% Real DeltaV
Delta_v_tot = 2e3; % [m/s]
DV_frac = 1; %%%%%%%%%%%%%%%%% DO NOT PLAY WITH THIS 

%%
H_min = max([DeltaV_Alt_Incl(DV,DV_frac) 2000])

%% Define Problem
% [A1.H, A1.theta0, A1.g_swath, A1.pho_x, A1.pho_gr, A1.num_phases]

X0 = [3000e3 deg2rad(20) 30e3 4 4 5];

lb = [H_min deg2rad(5) 5e3 2 2 5];
ub = [5000e3 deg2rad(65) 10e3 20 20 5];

nvars = length(X0);
pop_size = 100;
max_gen = nvars*200;

%rng default; % For Consistent results (same rand num seed at every execution)
opts = optimoptions(@paretosearch,'UseParallel', true, 'ParetoSetSize', pop_size,...
    'MaxIterations', max_gen, 'ConstraintTolerance', 1e-6,...% 'Display','iter',...
    'PlotFcn', {@psplotparetof},...% 'Display' , 'iter',...
    'ParetoSetChangeTolerance', 1e-6);    


%% Optimisation

[OptimParetoFront, OptimObjectiveVariables] = paretosearch(@ObjectiveFcnMultiobj,nvars,[],[],[],[],lb,ub,[],opts);
    
%% Refine Pareto Front

% % Filtering Pareto Objective Variables with Objective Variable Thresholds/Limits
% %[Final.HiRes_Cell Final.Aa Final.MissionTime A1.SNR_rx_dB]
% Thresh.RC = [1 8000]; % HiRes_Cell Voxel Thresholds [Min Max] [m^3]
% Thresh.Aa = [1 25]; % Antenna Area Thresholds [Min Max] [m^2]
% Thresh.MT = [1 15]; % Mission Time Thresholds [Min Max] [Yrs]
% 
% Front = [OptimObjectiveVariables(:,1) > Thresh.RC(1),...
%     OptimObjectiveVariables(:,2) > Thresh.Aa(1),...
%     OptimObjectiveVariables(:,3) > Thresh.MT(1)].*...
%     [OptimObjectiveVariables(:,1) < Thresh.RC(2),...
%     OptimObjectiveVariables(:,2) < Thresh.Aa(2),...
%     OptimObjectiveVariables(:,3) < Thresh.MT(2)];
% 
% Front2 = Front(:,1).*Front(:,2).*Front(:,3);
% Front3 = [Front2 Front2 Front2];
% FilteredFront = OptimObjectiveVariables.*Front3;

%% Show Optimum Condition
%[Final.HiRes_Cell Final.Aa Final.MissionTime]

figure
pointsize=50;
%scatter(FilteredFront(:,1), FilteredFront(:,2), pointsize, FilteredFront(:,3), 'filled');
scatter3(OptimObjectiveVariables(:,1), OptimObjectiveVariables(:,2), OptimObjectiveVariables(:,3), pointsize, OptimObjectiveVariables(:,4), 'filled');
title({'Filtered Optimum Objective Variables'; 'Corresponding to Optimum Pareto Front'})
xlim(gca,[0 10]) % Pixel Area
ylim(gca,[0 10]) % Antena Area
%set (gca, 'Color' , 'k' )
xlabel('Voxel Volume [m^3]')
ylabel('Antenna Area [m^2]')
zlabel('Mission Time [Years] ((GlobalScans HiRes = 2))')
cb = colorbar('eastoutside');
ylabel(cb, '-Rx SNR [-dB]' )

% %% Reverse solution
% [OptimCondition_idx, ~] = find(OptimObjectiveVariables(:,1) < voxel_limits(2) & OptimObjectiveVariables(:,1) > voxel_limits(1) & OptimObjectiveVariables(:,2) < Ant_a_limits(2) & OptimObjectiveVariables(:,2) > Ant_a_limits(1));
% format shortG
% OptimDesignVariables = OptimParetoFront(OptimCondition_idx, :);
% 
% % Optimum Antenna
% [F_ant.Lx, F_ant.Lz, F_ant.Pho_h, F_ant.Ptx, F_ant.delta_psi, F_ant.delta_theta] = Ant_design(OptimDesignVariables(4), OptimDesignVariables(5), OptimDesignVariables(3), OptimDesignVariables(1), OptimDesignVariables(2));
% 
% 
% %% Display Results
% disp(OptimDesignVariables)
% disp(F_ant)
% DisplayData(OptimDesignVariables);