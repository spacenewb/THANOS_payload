clc
clear 
close all


%%
%---------------------------SATURN ENVIROMENT-----------------------------%
T_ds=3;%[k]
alb_s=(0.34+0.49)*.5;%average btw bond and geometric abledo
alb_t=0.22; % titan albedo
F_a=1;% esempio
T_h=340.15;% hot case
T_c=250.15;% cold case
T_sc=300; %[K];


%---------------- SATURN RADIATIVE EMISSION---------------

T= 134;
% lamba=2898/T;%[micron]
% Wir=sigma*T^4
h=6.626e-34; %[j*s]
c=3e8; %[m/s]
K=1.38e-23; %[J/K]
f= @(x) (2.*pi.*h.*c.^2)./(x.^5).*(1./(exp(h.*c./(x.*T.*K))-1));

I_s=integral(f,7e-6,0.001);%Radiant emittance

%------------ TITAN RADIATIVE EMISSION--------------

 
T_t= 94;
% lamba=2898/T;%[micron]
sigma=5.67e-8;
% Wir=sigma*T^4
h=6.626e-34; %[j*s]
c=3e8; %[m/s]
K=1.38e-23; %[J/K]
f= @(x) (2.*pi.*h.*c.^2)./(x.^5).*(1./(exp(h.*c./(x.*T_t.*K))-1));

I_t=integral(f,7e-6,0.001);%Radiant emittance

%---------- VENUS EMITTANCE--------

T_ven= 698;
% lamba=2898/T;%[micron]
sigma=5.67e-8;
% Wir=sigma*T^4
h=6.626e-34; %[j*s]
c=3e8; %[m/s]
K=1.38e-23; %[J/K]
f= @(x) (2.*pi.*h.*c.^2)./(x.^5).*(1./(exp(h.*c./(x.*T_ven.*K))-1));

I_ven=integral(f,7e-6,0.001);%Radiant emittance


%%
%-----------------------------MATERIALS AND WALLS-------------------------%
Vol=5; % from baselines
F_ds=1;%view factor
L=2; % long size [m]
l=sqrt(Vol/L); %short size [m]
Ab=L*l;% large area of impact[m^2]
As=l*l;% small area of impact[m^2]
Asc_lat= L^2 - l^2;
Asc_fl= Ab;
Asc_ub=L^2;

% ALLUMINUM 2024 T-3
rho_al= 0.00277*10^6; %[Kg/m^3]
k_al=121.2;% conductivity[W/(m*K)] 
cp_al=880;%specific heat [J/(kg*K)]
alfa_al=0.24;%solar absorbance(visible)
epsilon_al=0.76;%thermal emittance(IR)
s_al=0.007;    %[m]
Gr_al=sigma*epsilon_al*F_ds;


% KAPTON 
rho_kap=0.00141*10^6; %[Kg/m^3]
k_kap=0.2;% conductivity[W/(m*K)]
cp_kap=1.423; %[J/kgK]
s_kap=0.002032/1000;   %[m] 0.08mil=0.002032
Gr_kap=sigma*0.83*F_ds;
% Ab_kap=L*l;% large area of impact[m^2]
% As_kap=l*l;% small area of impact[m^2]

% KAPTON ALLUMINIZED 0.08 MIL 
alfa_kap= 0.26;
epsilon_kap=0.24;

% DACRON DAC
rho_tfe= 1350; %[Kg/m^3]
k_tfe=0.2; % conductivity[W/(m*K)]


%--------------MLI PROVA---------
% l=1.58 metri

l_squagap=0.006;
n_squagap_v=180;
n_squagap_o=180;
n_squagap=n_squagap_v*n_squagap_o;
A_squagap=l_squagap^2; %[m^2]
A_net_new=As-n_squagap*A_squagap;

A_empty=As-A_net_new;

s_kap=0.002032/1000;
s_tfe=0.0008;

s_tot=10*s_tfe+11*s_kap

n_squagap_o_b=250;
n_squagap_b=n_squagap_v*n_squagap_o_b;
A_net_new_b=Ab-n_squagap_b*A_squagap;

A_empty_b=Ab-A_net_new_b;

Gr_kap_mli= sigma*((2-epsilon_kap)/epsilon_kap);
A_tot=4*Ab+2*As;
m_tot=4*(11*rho_kap*Ab*s_kap+10*rho_tfe*A_net_new_b*s_tfe)+2*(11*rho_kap*As*s_kap+10*rho_tfe*A_net_new*s_tfe);

m_kap=rho_kap*Ab*s_kap;

% MLI NET TFE
% l_squagap=0.005;
% n_squagap_v=90;
% n_squagap_o=180;
% n_squagap=n_squagap_v*n_squagap_o;
% A_squagap=l_squagap^2; %[m^2]
% A_pb=Ab;
% A_net=A_pb - (A_squagap*n_squagap);
% wire_th= (l-(n_squagap_v*l_squagap))/(n_squagap_v+1);

s_tfe_out=s_tfe; %[m]

% s_tfe_in=0.015/10; %[m]


 % THERMAL RESISTANCE [K/W]
% R_cond_al_b=s_al/(k_al*2*Ab);
% R_cond_al_s=s_al/(k_al*2*As);
% 
% R_cond_kap_b= s_kap/(k_kap*2*Ab);
% R_cond_kap_s= s_kap/(k_kap*2*As);
% 
% R_cond_tfe_b= s_tfe_out/(k_tfe*2*Ab);
% R_cond_tfe_s= s_tfe_out/(k_tfe*2*As);
% 
% R_rad_al_b=1/(sigma*epsilon_al*(Ab)*(T_h^2 + T_ds^2)*(T_h + T_ds));
% R_rad_al_s=1/(sigma*epsilon_al*(As)*(T_h^2 + T_ds^2)*(T_h + T_ds));
% 
% R_rad_kap_b=1/(sigma*epsilon_kap*(Ab)*(T_h^2 + T_ds^2)*(T_h + T_ds));
% R_rad_kap_s=1/(sigma*epsilon_kap*(As)*(T_h^2 + T_ds^2)*(T_h + T_ds));
% 
% R_rad_kap_bin=1/(sigma*epsilon_kap*(Ab)*(T_h^2 + T_sc^2)*(T_h + T_sc));
% R_rad_kap_sin=1/(sigma*epsilon_kap*(As)*(T_h^2 + T_sc^2)*(T_h + T_sc));
% 
% 
% R_mli_b=4*R_rad_kap_b + 3*R_cond_tfe_b + R_cond_al_b;
% R_mli_s=4*R_rad_kap_s + 3*R_cond_tfe_s + R_cond_al_s;
% 
% R_mli_bin=4*R_rad_kap_bin + 3*R_cond_tfe_b + R_cond_al_b;
% R_mli_sin=4*R_rad_kap_sin + 3*R_cond_tfe_s + R_cond_al_s;
% 
% R_cond_al_dis=s_al/(k_al*(As/2));
% R_rad_al_dis=1/(sigma*epsilon_al*(As/2)*(T_h^2 + T_ds^2)*(T_h + T_ds));
% 
% 
% R_PL_SC= (s_al/(2*Ab*k_al))+3*(1/(sigma*epsilon_kap*(2*Ab)*(T_h^2 + T_sc^2)*(T_h + T_sc)))+2*(s_tfe/(k_tfe*2*A_net));
% 
% % R_PL_SC= 2*(R_cond_al_b + 10*R_rad_kap_bin + 9*R_cond_tfe_b);
% R_PL_DS= 2*(R_cond_al_b + 10*R_rad_kap_b + 9*R_cond_tfe_b + R_cond_al_s + 10*R_rad_kap_s + 9*R_cond_tfe_s);

% MASSES

% v_al_b=Ab*s_al;
% v_kap_b=Ab*s_kap;
% v_tfe_out_b=Ab*s_tfe_out;
% v_tfe_in_b=Ab*s_tfe_in;
% 
% v_al_s=As*s_al;
% v_kap_s=As*s_kap;
% v_tfe_out_s=As*s_tfe_out;
% v_tfe_in_s=As*s_tfe_in;
% 
% m_al_b=rho_al*v_al_b;
% m_kap_b=rho_kap*v_kap_b;
% m_tfe_out_b=rho_tfe*v_tfe_out_b;
% m_tfe_in_b=rho_tfe*v_tfe_in_b;
% 
% m_al_s=rho_al*v_al_s;
% m_kap_s=rho_kap*v_kap_s;
% m_tfe_out_s=rho_tfe*v_tfe_out_s;
% 
% m_tot_al= 8*m_al_b + 4*m_al_s + 2*m_tfe_in_b + 2*m_tfe_out_b + 2*m_tfe_out_s;

m_rad=rho_al*s_al*As;


%------------------------HEAT INPUT COMPUTATION---------------------------% 

%-------- HOT CASE VENUS-----------
Wir_ven=I_ven;
Wsun_ven=2622;%[W/m^2]
alb_ven=(0.77+0.689)*0.5;

Qs_ven=Wsun_ven*alfa_al*Ab;
Qa_ven=Wsun_ven*alb_ven*Ab*alfa_al*F_a;
Qir_ven=Wir_ven*epsilon_al*Ab;

%------------- SATURN ENVIRONMENT----------

W_sun=14.82;% solar radiation[W/m^2]
W_ir_s=I_s;
W_ir_t=I_t;

% INSTRUMENT 
Qs_s=W_sun*alfa_al*Ab;%[W]
Qa_t=W_sun*alb_t*Ab*alfa_al*F_a;%[W]
Qir_t=W_ir_t*epsilon_al*Ab;%[W]
Qa_s=W_sun*alb_s*Ab*alfa_al*F_a;%[W]
Qir_s=W_ir_s*epsilon_al*Ab;%[W]
Q_ele=15; %[W]
Q_radar=60; %[W]
Q_in=Q_ele+Q_radar;
Q_ds=sigma*epsilon_al*F_ds*2*(Ab+As)*(T_h^4-T_ds^4);%[W]

% SPACE CRAFT
% Qs_s_sc=W_sun*alfa_al*Asc_fl;%[W]
% Qa_t_sc=W_sun*alb_t*Asc_fl*alfa_al*F_a;%[W]
% Qir_t_sc=W_ir_t*epsilon_al*Asc_fl;%[W]
% Qa_s_sc=W_sun*alb_s*Asc_fl*alfa_al*F_a;%[W]
% Qir_s_sc=W_ir_s*epsilon_al*Asc_fl;%[W]








%% 
%---------------------------ECLIPSES CONFIGURATION------------------------%
r=1500; %orbit
R_t=2576; %titan
R_s=58232; %saturn
D_st=1224526; %distance Saturn-Titan(circular)

D_sct=r+R_t;
beta=asind(R_t/D_sct);
D_scs=r+R_t+D_st;
gamma=asind(R_s/D_scs);

%%

Sat_Surf= 4.2*10^16; %[m^2] Surface of Saturn
diam_s=120536; %[Km]
R_s=diam_s/2;  %[Km]
R_t=2576; %[Km]
t_simul=365*24*3600; %[s]


% SATURN ORBIT

mu_sun=132712440018; %[Km^3/s^2]
peri_s=1352550000;    %[Km]
apo_s=1515500000;      %[Km]
semajax_s=(apo_s+peri_s)/2;  %[Km]
ecc_s= (apo_s-peri_s)/(apo_s + peri_s); %[/]
period_s=29.45*365*24*3600; %[s]
mean_motion_s= (2*pi)/period_s; %[rad/s]
incl_s=deg2rad(2.485); %[rad]
incl_ax_s=deg2rad(26.73); % [rad]
mean_r_s=1.427*10^9; %[Km]



% TITAN ORBIT
mu_s=37931187; %[Km^3/s^2]
semajax_t=1221830; %[Km]
ecc_t=0.0293; %[/]
period_t=15.945421*24*3600; %[s]
incl_t=deg2rad(0.3485); %[rad]
mean_motion_t= (2*pi)/period_t; %[rad/s]
incl_t_ecl=incl_t+incl_s;

% S/C ORBIT
mu_t=8971.15; %[Km^3/s^2]
semajax_sc=1500+R_t; %[Km]
ecc_sc=0; %[/]
period_sc= 2*pi*sqrt(semajax_sc^3/mu_t); %[s]
incl_sc=deg2rad(90+45); %[]
mean_motion_sc= 2*pi/period_sc;

%max_step is the integration max step to input in simulink when the 
%integrator is characterized. It is important to set this value small 
%enough wrt the integration time otherwhise the eclipse block could miss 
%some points of integration, thus some eclipses are not computed. 
max_step=5000;

%-------------- ANTENNA+REFLECTOR--------

l_ant=1.31;%[m]
w_ant=0.216;

l_ref=5;
d_ref=3;
l_ap_ref=5;
w_ap_ref=0.28086; %aperture width
w_cl_ref=0.625; %closed width
t_ref=0.625; % thikness

m_ant=0.964367; % antenna mass
m_ref=50.34085; % reflector mass

A_ref=pi*l_ref*d_ref/4;
A_ant=l_ant*w_ant;

F_ar=0.8;
F_ra=F_ar*A_ant/A_ref;

Gr_ar=sigma/((1-epsilon_al)/(epsilon_al*A_ant)+1/(F_ar*A_ant)+(1-epsilon_al)/(epsilon_al*A_ref));
Gr_ra=sigma/((1-epsilon_al)/(epsilon_al*A_ref)+1/(F_ra*A_ref)+(1-epsilon_al)/(epsilon_al*A_ant));

d_struc=0.0762;
A_struc=pi/4*d_struc^2;
