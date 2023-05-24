clear all
clc

H = 900e3;
total_swath = 30e3;
sat_bend_angle = deg2rad(10)

pho_x = 2.5;
pho_gr = 2.5;

b = H*tan(sat_bend_angle);

column_elements = 5;
n = 1:column_elements;

equiangle_beam_angle = ones(1, column_elements).*rad2deg((atan((total_swath+H*tan(sat_bend_angle))/H)-sat_bend_angle)/column_elements) % Resolution changes?
equiswath_beam_angles = rad2deg(atan((H*tan(sat_bend_angle) + n.*(total_swath/column_elements))/H)-sat_bend_angle);


for i=2:column_elements
    equiswath_beam_angles(i) = equiswath_beam_angles(i) - sum(equiswath_beam_angles(1:(i-1)));
    ant_look_angles(i) = sat_bend_angle + sum(equiswath_beam_angles(1:(i-1))) + equiswath_beam_angles(i)/2;
end
equiswath_beam_angles
ant_look_angles(1) = sat_bend_angle + equiswath_beam_angles(1)/2;

for i=1:column_elements
    P_arr(i) = array_calc(deg2rad(ant_look_angles(i)),pho_x,pho_gr,deg2rad(equiswath_beam_angles(i)));
end

P_arr
P_total = sum(P_arr)