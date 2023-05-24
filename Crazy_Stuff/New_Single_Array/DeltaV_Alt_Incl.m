function [H_T] = DeltaV_Alt_Incl(Delta_v_tot,DV_frac)

Delta_v_tot = Delta_v_tot/1000; % Conversion from m/s to km/s

%Planet Feature
R_T=2574.7;                         %[Km]
M_T=1.34*(10^23);                   %[Kg]
mu_T= M_T * 6.67*(10^-20);          %[Km^3/s^2]

%Plane change man. estimation
%Based on the total Delta v given as constraint, number of plane change
%maneouvres desired and the total angle to spann the computation of step of
%inclination, altitude of the orbit and propellant mass are done.
%Assumption based on Cassini are taken for Isp and total mass.

%input for maneouvre
Delta_v_tot=Delta_v_tot*DV_frac;    %[Km/s]
n_man=2;                            %[45-90/90-135]
obj_inc=90;                         %[deg]

%values
step_inc=obj_inc/n_man; 
delta=deg2rad(step_inc);            %[rad]
Delta_v=Delta_v_tot/n_man ;          %[Km/s]
v_circ=Delta_v/(2*sin(delta/2));    %[Km/s]
r_T = mu_T/v_circ^2;                %[Km]
a_T=r_T;                            %[Km] 
H_T = a_T - R_T;                      %[Km]
H_T = H_T*1000;                     % [m]
end
