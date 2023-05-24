function [ObjectiveVariables] = EvalFcn(DesignVariables)
%% Input Design Variables
OrbAlt = DesignVariables(1);
Theta0 = DesignVariables(2);
g_swath = DesignVariables(3);
Res_Ant1.R = DesignVariables(4);
Res_Ant1.Az = DesignVariables(5);
Res_Ant2.R = DesignVariables(6);
Res_Ant2.Az = DesignVariables(7);

% OrbAlt = 900e3; % [m]
% Theta0 = 25; % [Rad]
% g_swath = 15e3; % Manual Ground Swath
% Res_Ant1.R = 1; % [m]
% Res_Ant1.Az = 4; % [m]
% Res_Ant2.R = 5; % [m]
% Res_Ant2.Az = 2; % [m]

%% Design Constraints
P_gen = 110; % Power Source in [W] (RTG)

g_swath_manual = 1; % Switch to select ground swath [0: Automatic, 1: Manual]
g_swath_ratio = 4; % Ratio of g_swath of Ant1/Ant2

GlobalScans_lo = 2; % Min number of LowRes global scans during mission (Maybe Not Required)
GlobalScans_hi = 2; % Min number of HiRes global scans during mission (Maybe Sufficient)

%% Orbital
% Circular Polar Orbit Assumption
% Ground Swath is Defined by Orbital Alt (deltaLambda), Change switch if required!!!!!!

Orbital.G = 6.67408e-11; % Standard Gravitational Constant
Orbital.Eyr = 365.24; % Earth Year in Earth Days [days]
Orbital.Eday = 86400.002; % Earth Day in [s]

Orbital.Rb = 2.57473e6; % Radius of Titan
Orbital.Mb = 1.3452e23; % Mass of Titan
Orbital.Wb = 4.5608e-06; % Angular Rotation Speed of Titan [Rad/sec]
Orbital.MU = Orbital.G*Orbital.Mb; % Gravitational Parameter of Titan

Orbital.H = OrbAlt; % Orbital Altitude [m]
Orbital.a = Orbital.H + Orbital.Rb; % Orbital SMA [m]
Orbital.Vel = (Orbital.MU/Orbital.a)^0.5; % Orbital velocity [m/s]
Orbital.T = (2*pi*(((Orbital.Rb+Orbital.H)^3)/Orbital.MU)^0.5); % Circular Orbit [s]
Orbital.dLambda = Orbital.Wb*Orbital.T; % Delta Lambda of Ground Track [Rad]
Orbital.dLSwath = Orbital.dLambda*Orbital.Rb; % GroundSwath of DeltaLambda

%% Antenna
% Ant1: High Range Resolution
% Ant2: High Azimuth Resolution

% Input
Ant1.theta0 = Theta0; % Elevation Angle of Ant1
Ant2.theta0 = Theta0; % Elevation Angle of Ant2

Ant1.res = [Res_Ant1.R Res_Ant1.Az]; % Range & Azimuth Resolution of Ant1
Ant2.res = [Res_Ant2.R Res_Ant2.Az]; % Range & Azimuth Resolution of Ant2

if g_swath_manual > 0
    Ant1.g_swath = g_swath; % Ground Swath of Ant1 (Manually Set)
else
    Ant1.g_swath = Orbital.dLSwath; % Ground Swath of Ant1 (Based on Orbital.dLSwath)
end

Ant2.g_swath = Ant1.g_swath/g_swath_ratio; % Ground Swath of Ant2 (< Ant1)

% Calculated
[Ant1.Lx, Ant1.Lz, Ant1.pho_h, Ant1.Ptx] = Ant_design(Ant1.res(2), Ant1.res(1), Ant1.g_swath, Orbital.H, Ant1.theta0);
[Ant2.Lx, Ant2.Lz, Ant2.pho_h, Ant2.Ptx] = Ant_design(Ant2.res(2), Ant2.res(1), Ant2.g_swath, Orbital.H, Ant2.theta0);

Ant1.Aa = Ant1.Lx*Ant1.Lz; % Antenna Area of Ant1
Ant2.Aa = Ant2.Lx*Ant2.Lz; % Antenna Area of Ant1

%% Energy Budget
E_gen_orb = P_gen*Orbital.T; % Energy Generated in 1 Orbit

%% Scan Parameters
Scan.res_low = Ant1.res; % Range and Azimuth of Ant1
Scan.res_hi = [Ant1.res(1) Ant2.res(2)]; % Range and Azimuth of combined best

Scan.ToS_orb = E_gen_orb/(Ant1.Ptx+Ant2.Ptx); % Maximum Scanning Time in 1 Orbit [s]
Scan.DC_orb = Scan.ToS_orb/Orbital.T; % Duty Cycle of SAR in 1 Orbit
Scan.L_Az_orb = Scan.ToS_orb*Orbital.Vel; % Length of Scanning in each Orbit along Azimuth [m]
Scan.A_orb_lo = Scan.L_Az_orb*Ant1.g_swath; % Area of Surface Scanned in LowRes every Orbit [m^2]
Scan.A_orb_hi = Scan.L_Az_orb*Ant2.g_swath; % Area of Surface Scanned in HiRes every Orbit [m^2]

if g_swath_manual > 0
    Scan.Revisit_N = (2*pi/Orbital.dLambda)*(Orbital.dLSwath/Ant1.g_swath); % No. of Orbits to revisit same ground track
else
    Scan.Revisit_N = (2*pi/Orbital.dLambda); % No. of Orbits to revisit same ground track
end

Scan.Revisit_T = Scan.Revisit_N*Orbital.T; % Time for revisiting same ground track

Scan.Revisits = ceil(1/Scan.DC_orb); % Approx. No. of Revisits req. to cover same ground track

Scan.Global_T_lo = Scan.Revisits*Scan.Revisit_T; % Time to globally scan Titan Once LowRes
Scan.Global_T_hi = Scan.Global_T_lo*g_swath_ratio; % Time to globally scan Titan Once HiRes

Scan.Mission_lo = max(GlobalScans_hi*g_swath_ratio, GlobalScans_lo); % Total No. of LowRes Global Scans in Mission
Scan.Mission_T = Scan.Mission_lo*Scan.Global_T_lo; % Total Mission Time [s]

%% Output

% Computations
Final.Ptx = Ant1.Ptx+Ant2.Ptx; % [w]
Final.Aa = Ant1.Aa+Ant2.Aa; % [m^2]
Final.HiRes_Cell = Scan.res_hi(1)*Scan.res_hi(2); % [m^2]
Final.MissionTime = (Scan.Mission_T)/Orbital.Eday/Orbital.Eyr; % In Earth Years

% Objective Variable
ObjectiveVariables = [Final.HiRes_Cell Final.Aa Final.MissionTime];

end
