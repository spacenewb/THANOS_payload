function[c, ceq] = ConstraintFcn(X)
%Ant1.ResR>Ant2.ResR
%X(6)-X(4)<0

%Ant1.ResR>Ant1.ResAz
%X(5)-X(4)<0

%Ant2.ResAz>Ant2.ResR
%X(6)-X(7)<0

%%
c = [X(6)-X(4); X(5)-X(4); X(6)-X(7)];

%constraint on pixel aspect ratio for each antenna
res1_AR = 2;
res2_AR = 3;

ceq = [(X(4)/X(5))-res1_AR; (X(7)/X(6))-res2_AR];

%%
% c = [X(6)-X(4); X(5)-X(4); X(6)-X(7)];
% 
% ceq = [];

end
