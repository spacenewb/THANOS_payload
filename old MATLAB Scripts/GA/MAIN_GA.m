clc
clear all

%% Nomenclature
%DesignVariables: The variables that we use to design the SAR eg:Alt,swath,elevation..
%ObjectiveVariable: The variable used to judge the optimisation score eg:power,mass..
%ParetoFront: The set of Design variables tht are a result of Pareto optimum ObjectiveVariables

%% Switch to choose Pareto Opt or Single Variable Opt
MultiObjSwitch = 0; % Use non-zero to perform MultiObjective Pareto GA Optimisation

%% Create Parallel Pool
if gcp('nocreate')==0
parpool()
end

%% Define Problem
%(OrbAlt, Theta0, g_swath, Res_Ant1.R, Res_Ant1.Az, Res_Ant2.R, Res_Ant2.Az)

X0 = [900e3 deg2rad(20) 30e3 1 4 3 2];

lb = [900e3 deg2rad(0) 10e3 10 10 10 10];
ub = [1500e3 deg2rad(35) 50e3 100 100 100 100];

nvars = length(X0);

%rng default; % For Consistent results (same rand num seed at every execution)

if MultiObjSwitch == 0
    opts = optimoptions(@ga,'UseParallel', true, 'PopulationSize', 500,...
        'MaxGenerations', nvars*200, 'MaxStallGenerations', nvars*20,...
        'ConstraintTolerance', 1e-6, 'CrossoverFraction', 0.8,...
        'PlotFcn', {@gaplotbestf,@gaplotscorediversity}, 'Display' , 'iter',...
        'FunctionTolerance', 1e-6);
else
    opts = optimoptions(@gamultiobj,'UseParallel', true, 'PopulationSize', 500,...
        'MaxGenerations', nvars*200, 'MaxStallGenerations', nvars*20,...
        'ConstraintTolerance', 1e-6, 'CrossoverFraction', 0.8,...
        'PlotFcn', {@gaplotscorediversity},...% 'Display' , 'iter',...
        'FunctionTolerance', 1e-6);    
end


%% Optimisation
if MultiObjSwitch == 0
    [OptimDesignVariables, OptimObjectiveVariable] = ga(@ObjectiveFcn,nvars,[],[],[],[],lb,ub,[],opts);
else
    [OptimParetoFront, OptimObjectiveVariables] = gamultiobj(@ObjectiveFcnMultiobj,nvars,[],[],[],[],lb,ub,[],opts);
end

%% Refine Pareto Front
if MultiObjSwitch > 0
    % Filtering Pareto Objective Variables with Objective Variable Thresholds/Limits
    %[Final.HiRes_Cell Final.Aa Final.MissionTime]
    Thresh.RC = [1 400]; % HiRes_Cell Thresholds [Min Max] [m^2]
    Thresh.Aa = [1 15]; % Antenna Area Thresholds [Min Max] [m^2]
    Thresh.MT = [4 10]; % Mission Time Thresholds [Min Max] [Yrs]

%     Thresh.Max = [400 15 10];
%     Thresh.Min = [1 1 4];
% 
%     NumObjVars = length(OptimObjectiveVariables(1,:));
%     filter = zeros(size(OptimObjectiveVariables));
%     for i=1:NumObjVars
%         filter(:,i) = OptimObjectiveVariables(:,i)>=Thresh.Min(i) & OptimObjectiveVariables(:,i)<=Thresh.Max(i);
%         if i == NumObjVars
%             Filter = filter(:,1);
%             for j=1:(NumObjVars-1)
%                 Filter = Filter.*filter(:,j+1);
%             end
%             FFilter = repmat(Filter, 1, NumObjVars);
%         end
%     end
%     FilteredFront = OptimObjectiveVariables.*FFilter;


    Front = [OptimObjectiveVariables(:,1) > Thresh.RC(1),...
        OptimObjectiveVariables(:,2) > Thresh.Aa(1),...
        OptimObjectiveVariables(:,3) > Thresh.MT(1)].*...
        [OptimObjectiveVariables(:,1) < Thresh.RC(2),...
        OptimObjectiveVariables(:,2) < Thresh.Aa(2),...
        OptimObjectiveVariables(:,3) < Thresh.MT(2)];
    
    Front2 = Front(:,1).*Front(:,2).*Front(:,3);
    Front3 = [Front2 Front2 Front2];
    FilteredFront = OptimObjectiveVariables.*Front3;
end

%% Show Optimum Condition
%[Final.HiRes_Cell Final.Aa Final.MissionTime]
if MultiObjSwitch == 0
    OptimumSoln = EvalFcn(OptimDesignVariables);
    disp(OptimDesignVariables)
    disp(OptimumSoln)
else
    figure(1)
    pointsize=50;
    scatter(FilteredFront(:,1), FilteredFront(:,2), pointsize, FilteredFront(:,3), 'filled');
    title({'Filtered Optimum Objective Variables'; 'Corresponding to Optimum Pareto Front'})
    xlim(gca,[0 10]) % Pixel Area
    ylim(gca,[0 10]) % Antena Area
    set (gca, 'Color' , 'k' )
    xlabel('HiRes Cell Area [m^2]')
    ylabel('Antenna Area [m^2]')
    cb = colorbar('eastoutside');
    ylabel(cb, 'Mission Time [Years] ((GlobalScans HiRes = 2))' )
end

%% Optimum Antenna
[HRR.Lx, HRR.Lz, HRR.PhoX, HRR.Ptx] = Ant_design(OptimDesignVariables(5), OptimDesignVariables(4), OptimDesignVariables(3), OptimDesignVariables(1), OptimDesignVariables(2));
[HAR.Lx, HAR.Lz, HAR.PhoX, HAR.Ptx] = Ant_design(OptimDesignVariables(7), OptimDesignVariables(6), OptimDesignVariables(3), OptimDesignVariables(1), OptimDesignVariables(2));
%% figure

% figure(2)
% pointsize=50;
% scatter(OptimObjectiveVariables(:,3), OptimObjectiveVariables(:,2), pointsize, OptimObjectiveVariables(:,1), 'filled');
% title({'Filtered Optimum Objective Variables'; 'Corresponding to Optimum Pareto Front'})
% xlim(gca,[0 10]) % Mission Time
% ylim(gca,[0 10]) % Antenna Area
% set (gca, 'Color' , 'k' )
% xlabel('Mission Time [Years] ((GlobalScans HiRes = 2))' )
% ylabel('Antenna Area [m^2]')
% cb = colorbar('eastoutside');
% ylabel(cb, 'HiRes Cell Area [m^2]')

%% fig

%figure(3)
%plot(OptimObjectiveVariables(:,2),OptimObjectiveVariables(:,3),OptimObjectiveVariables(:,2),OptimObjectiveVariables(:,1),'Marker')