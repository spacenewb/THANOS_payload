clc
clear all
close all

%% PHYSICAL CONSTANTS
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant

%%


%% For the range antenna: A1 (Range Res High)

A1.H = 1500e3;
A1.theta0 = deg2rad(10);
A1.g_swath = 30e3;
A1.pho_x = 5;
A1.pho_gr = 2.5;

[A1.Lx, A1.Lz, A1.pho_h, A1.Ptx, A1.delta_psi, A1.delta_theta] = Ant_design(A1.pho_x, A1.pho_gr, A1.g_swath, A1.H, A1.theta0);

%% For the azimuth antenna: A2 (Az Res High)

A2.H = 900e3;
A2.theta0 = 0.0000000001; % Only scanning at boresight along the swath azimuth? Or sat change attitude?
A2.Lx = A1.Lz;
A2.Lz = A1.Lx;

A2.pho_gr = A1.pho_x; % Can change this 
A2.pho_gr = 3; % Can change this 

A2.g_swath = ((A2.H/cos(A2.theta0))*A1.delta_psi)/cos(A2.theta0); % The new ground swath depends on the Azimuth Beam width of prev. Config.

[A2.pho_x, A2.pho_gr, A2.pho_h, A2.Ptx, A2.delta_psi, A2.delta_theta] = Ant_design_rev(A2.Lx, A2.Lz, A2.pho_gr, A2.g_swath, A2.H, A2.theta0);


%%

Total_pow = A1.Ptx + A2.Ptx;
GS_ratio = A1.g_swath/A2.g_swath;

A1.NO = 2*pi*R_B/A1.g_swath;
A2.NO = 2*pi*R_B/A2.g_swath;


A1.T_orbit = (2*pi*(((R_B+A1.H)^3)/Mu_B)^0.5)/86400; % Earth Days
A2.T_orbit = (2*pi*(((R_B+A1.H)^3)/Mu_B)^0.5)/86400; % Earth Days

power_budget = 110; % Baseline [W]

A1.DC = ((power_budget*A1.T_orbit*86400)/A1.Ptx)/(A1.T_orbit*86400);
A2.DC = ((power_budget*A2.T_orbit*86400)/A2.Ptx)/(A2.T_orbit*86400);

A1.MD = ceil(1/A1.DC)*A1.NO*A1.T_orbit; % In Days
A2.MD = ceil(1/A2.DC)*A2.NO*A2.T_orbit;

Res.GR = [A1.pho_gr A2.pho_gr];
Res.Az = [A1.pho_x A2.pho_x];

Res.Hi1 = [A1.pho_gr A2.pho_x];
Res.Hi2 = [A2.pho_gr A1.pho_x];

%%

[Res_Cell_A, Res.idx] = min([(Res.Hi1(1)*Res.Hi1(2)) (Res.Hi2(1)*Res.Hi2(2))]);

if Res.idx ==1
    Res.Hi = Res.Hi1;
else
    Res.Hi = Res.Hi2;
end

MD = (A1.MD + A2.MD)/365.25 % Years
Ant_A = A1.Lx*A1.Lz; %
Ant_Dim = [A1.Lx A1.Lz]
Hi_Resoluton_Cell = Res.Hi
Res_Cell_A
Res_Alt = [A1.pho_h A2.pho_h]
MD_Individual = [A1.MD A2.MD]./365.25
Ptx_Power = [A1.Ptx A2.Ptx]

Long = linspace(0,359,360);
for i=1:360
    Decorr_Time(i) = ((Long(i)/360)*((A2.MD-A1.MD)/365.25))+(A1.MD/365.25); % Years
end

plot(Long, Decorr_Time)

%% Display
% A1
% A2
