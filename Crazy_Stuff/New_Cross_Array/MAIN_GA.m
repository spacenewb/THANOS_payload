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

%% Define Problem
%(OrbAlt, Theta0, g_swath, Res_Ant1.R, Res_Ant1.Az, Res_Ant2.R, Res_Ant2.Az)

X0 = [900e3 deg2rad(20) 30e3 1 4 3 2];

lb = [900e3 deg2rad(0) 10e3 2 3 3 2];
ub = [3000e3 deg2rad(35) 50e3 20 30 30 20];

nvars = length(X0);

%rng default; % For Consistent results (same rand num seed at every execution)

opts = optimoptions(@gamultiobj,'UseParallel', true, 'PopulationSize', 500,...
    'MaxGenerations', nvars*200, 'MaxStallGenerations', nvars*20,...
    'ConstraintTolerance', 1e-6, 'CrossoverFraction', 0.8,...
    'PlotFcn', {@gaplotscorediversity},...% 'Display' , 'iter',...
    'FunctionTolerance', 1e-6);    

%% Optimisation
[OptimParetoFront, OptimObjectiveVariables] = gamultiobj(@ObjectiveFcnMultiobj,nvars,[],[],[],[],lb,ub,[],opts);

% %% Refine Pareto Front
% % Filtering Pareto Objective Variables with Objective Variable Thresholds/Limits
% %[Final.HiRes_Cell Final.Aa Final.MissionTime]
% Thresh.RC = [1 8000]; % HiRes_Cell Voxel Thresholds [Min Max] [m^3]
% Thresh.Aa = [1 25]; % Antenna Area Thresholds [Min Max] [m^2]
% Thresh.MT = [3 15]; % Mission Time Thresholds [Min Max] [Yrs]
% 
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

figure(1)
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
% %% Optimum Antenna
% %[F_ant.Lx, F_ant.Lz, F_ant.Pho_h, F_ant.Ptx, F_ant.delta_psi, F_ant.delta_theta] = Ant_design(OptimDesignVariables(4), OptimDesignVariables(5), OptimDesignVariables(3), OptimDesignVariables(1), OptimDesignVariables(2));
% % Optimum Antenna
% [HRR.Lx, HRR.Lz, HRR.PhoX, HRR.Ptx, HRR.delta_psi, HRR.delta_theta] = Ant_design(OptimDesignVariables(5), OptimDesignVariables(4), OptimDesignVariables(3), OptimDesignVariables(1), OptimDesignVariables(2));
% [HAR.Lx, HAR.Lz, HAR.PhoX, HAR.Ptx, HAR.delta_psi, HAR.delta_theta] = Ant_design(OptimDesignVariables(7), OptimDesignVariables(6), OptimDesignVariables(3), OptimDesignVariables(1), OptimDesignVariables(2));
% 
% %% Display Results
% disp(OptimDesignVariables)
% disp(HRR)
% disp(HAR)
% DisplayData(OptimDesignVariables);