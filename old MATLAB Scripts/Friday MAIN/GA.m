clc
clear all

%(pho_x, pho_gr, g_swath, H, theta0)
X0 = [10 10 50e3 900e3 pi/8];
lb = [1 1 10e3 900e3 deg2rad(5)];
ub = [100 100 50e3 1500e3 deg2rad(45)];

options = optimoptions ('ga', 'PlotFcn' , {@ gaplotbestf}, ... 
    'Display' , 'iter' );
options.InitialPopulationMatrix = X0;

[x, fval] = ga(@ObjectiveFcn,5,[],[],[],[],lb,ub,[],options);

x
fval;

angle = rad2deg(x(5));
[Lx, Lz, pho_h, G_ant_db, Ptx, N_coverage] = Ant_design2(x(1),x(2),x(3),x(4),x(5));