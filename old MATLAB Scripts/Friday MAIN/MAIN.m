clear all
clc

%%%%%% Only edit the parameters specified in "Assumptions" section %%%%%%

%% INITIALISATION
global body_mass Grav_const c R_B Mu_B K
global f0 sigma0_NESZ F sigma_0 eta duty_cycle
global H v_orb theta0

%% PHYSICAL CONSTANTS
body_mass = 1.3452e23; % Mass of the parent body
Grav_const = 6.67408e-11; % Standard Gravitational Constant
c = physconst('LightSpeed'); % Wave propagation velocity [m/s]
R_B = 2.57473e6; % Radius of Body (Earth = physconst('EarthRadius'))
Mu_B = Grav_const*body_mass; % Gravitational Parameter of Body (3.986e14 for earth)
K = physconst('boltzmann'); % Boltzmann constant

%% GIVEN DATA
f0 = 13.78e9; % Carrier frequency [Hz]
sigma0_NESZ = 10^(0.1*(-22)); % Noise Equivalent Sigma Zero [dB]
F = 10^(0.1*(5)); % Reciever Noise Figure [dB]
sigma_0 = 10^(0.1*(-15)); % Sigma Zero [dB]
eta = 0.7; % Antenna Efficiency
duty_cycle = 0.3; % Duty Cycle

%% ASSUMPTIONS
H = 1500e3; % orbit height [m]
v_orb = (Mu_B/(R_B+H))^0.5; % orbit velocity [m/s]
theta0 = deg2rad(5); % Elevation Angle [rad]

res.AR = 2; % Low Resolution/ High Resolution
res.low = 5; % Low Resolution
res.hi = res.low/res.AR; % High Resolution

swath_ratio = 5; %HRR/HAR
HRR.g_swath = 47e3;
HAR.g_swath = HRR.g_swath/swath_ratio;

power_budget = 110; % Baseline [W]

%% ANTENNAE COMPUTATIONS
% HRR: High Range Resolution Antenna
% HAR: High Azimuth Resolution Antenna

[HRR.Lx, HRR.Lz, HRR.pho_h, HRR.G_ant_db, HRR.Ptx] = Ant_design(res.low, res.hi, HRR.g_swath);
[HAR.Lx, HAR.Lz, HAR.pho_h, HAR.G_ant_db, HAR.Ptx] = Ant_design(res.hi, res.low, HAR.g_swath);

%% ORBITAL CALC
% Circular 0 Inclination Orbit Assumption for Orbital Time Calculation
% No Downtime considered for downlink, power, orbit change etc.
% Sun-synchrnous not considered

HRR.N_coverage = pi/asin(HRR.g_swath/(2*R_B)); % Number of orbits to cover whole surface

HRR.T_orbit = (2*pi*(((R_B+H)^3)/Mu_B)^0.5)/86400; % Earth Days
HRR.T_revisit = HRR.N_coverage*HRR.T_orbit; % Earth Days
HAR.T_revisit = HRR.T_revisit*swath_ratio; % Earth Days


%% OUTPUT
HRR;
HAR;

Alt_resolution = ["HRR" ,HRR.pho_h; "HAR" ,HAR.pho_h] % Altitude Resolution of Both Antennae [m]

Ant_footprint = ["Combined" , max(HRR.Lx,HAR.Lx), max(HRR.Lz,HAR.Lz); ...
    "HRR" , HRR.Lx, HRR.Lz; "HAR", HAR.Lx, HAR.Lz] % Antenna Footprint [m] of Both Antennae & Combined

P_total = sum([HRR.Ptx, HAR.Ptx]) % Combined Antenna Power [W]

Mission_Duration = [HRR.T_revisit, HAR.T_revisit] % Earth Days for Low-Res and High-Res respectively

%% Optimisation

GS_vec = linspace(5e3,70e3,10001);

for i=1:length(GS_vec)
    
    g_swath1 = GS_vec(i);
    g_swath2 = g_swath1/swath_ratio;
    
    [~,~,~,~,p1] = Ant_design(res.low, res.hi, g_swath1);
    [~,~,~,~,p2] = Ant_design(res.hi, res.low, g_swath2);
    
    N_cov(i) = pi/asin(g_swath1/(2*R_B)); % Number of orbits to cover whole surface Mid Resolution
    
    P_vec(i) = p1 + p2;
    
    energy_budget_J = (power_budget*HRR.T_orbit*86400);
    
    ToA = (energy_budget_J)/P_vec(i);
    dut_cyc(i) = ToA/(HRR.T_orbit*86400);
    
    n_cov(i) = ceil(1/dut_cyc(i))*N_cov(i)*HRR.T_orbit;
    
    areaOfScan(i) = ToA*v_orb*g_swath1;
    
end
%%
figure
y = linspace(0,4000);
x = 50.*ones(1,length(y));
plot(GS_vec/1000, P_vec, GS_vec/1000, n_cov, GS_vec/1000, areaOfScan/1e9, GS_vec/1000, dut_cyc*1e2, x, y, '-.k','LineWidth', 2)
xlabel('Ground Swath [Km]')
ylabel('Magnitude')
legend(["Antenna Power [W]" "Global Scan Time [Days]" "Area of 1 scan [1e-3 Km^2]" "Duty Cycle [1e2]"])
