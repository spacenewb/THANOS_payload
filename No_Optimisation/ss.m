function [Ant_Mass] = Ant_Mass_eval(Acq_Modes, Ant_Area, Ant_res)
data1 = load('MassAreaDensity.mat');
Modes1 = table2array(data1.MassAreaDensity(:,1));
MassDensity1 = table2array(data1.MassAreaDensity(:,2));
MassD1 = interp1(Modes1, MassDensity1, Acq_Modes, 'makima')
Mass1 = MassD1*Ant_Area

data2 = load('Mass_Modes.mat');
Modes2 = table2array(data2.MassModes(:,1));
Massvals2 = table2array(data2.MassModes(:,2));
Mass2 = interp1(Modes2, Massvals2, Acq_Modes, 'makima')

data3 = load('ResMass.mat');
Ant_res = Ant_res;
Masses3 = table2array(data3.ResMass(:,1));
Resolution = table2array(data3.ResMass(:,2));
Mass3 = interp1(Resolution, Masses3, Ant_res, 'linear')

masses = [Mass1 Mass2 Mass3];
AntMass = sum(masses)/numel(masses)
end