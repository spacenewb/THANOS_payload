function[ObjectiveVariables] = ObjectiveFcnMultiobj(DesignVariables)

A1.H = DesignVariables(1);
A1.theta0 = DesignVariables(2);
A1.g_swath = DesignVariables(3);
A1.pho_x = DesignVariables(4); 
A1.pho_gr = DesignVariables(5);
A1.num_phases = DesignVariables(6);

%% PHYSICAL CONSTANTS
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant

power_budget = 110; % Baseline [W]

%% Design Antenna Phase A1, A2, C1, C2 & Polar

% A1.H = H_min;
% A1.theta0 = deg2rad(30);
% A1.g_swath = 30e3;
% A1.pho_x = 5; 
% A1.pho_gr = 2;

HiRes = min([A1.pho_x A1.pho_gr]);

[A1.Lx, A1.Lz, A1.pho_h, A1.Ptx, A1.delta_psi, A1.delta_theta, A1.SNR_rx_dB] = Ant_design(A1.pho_x, A1.pho_gr, A1.g_swath, A1.H, A1.theta0);

% Compute Mission Parameters
A1.NO = 2*pi*R_B/A1.g_swath;
A1.T_orbit = (2*pi*(((R_B+A1.H)^3)/Mu_B)^0.5)/86400; % Earth Days
A1.DC = ((power_budget*A1.T_orbit*86400)/A1.Ptx)/(A1.T_orbit*86400);
A1.MD = ceil(1/A1.DC)*A1.NO*A1.T_orbit; % In Days

%% Final Mission Parameters
Final.NO = A1.num_phases*A1.NO;
Final.MD = A1.num_phases*A1.MD/365.25; % Years
Final.HiRes = [HiRes HiRes A1.pho_h];
Final.ResVol = Final.HiRes(1)*Final.HiRes(2)*Final.HiRes(3);
Final.AntArea = A1.Lx*A1.Lz;

ObjectiveVariables = [Final.ResVol Final.AntArea Final.MD -A1.SNR_rx_dB]; % Final.Error];

end
