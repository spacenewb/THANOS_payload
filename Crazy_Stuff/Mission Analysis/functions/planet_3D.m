function planet_3D(radius,dist,image_in)
% Plot Planet 3-D
%
% Inputs:
% image_in      : image of Titan to plot over the sphere
%
% Plots:
% Planet 3-D
% 
% Contributors:
% Gaballo Paolo
%  
% Versions:
% 2021-04-13, third version
% 2021-04-12, second version (Titan_3D)
% 2021-02-13, first version (earth)


if nargin == 2
    image_in = 'Titan_plan.jpg';
end

Rt = radius;  

% Define the number of panels to be used to model the sphere 
npanels = 180;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(dist, 0, 0, Rt, Rt, Rt, npanels);

% Create the globe with the surf function
globe = surf(x, y, -z, 'EdgeColor', 'none', 'FaceColor', 'none');
hold on;
axis equal;
%Earth_image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

cdata = imread(image_in);
% cdata = imread('EarthTexture.jpg');

set(globe,'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha',1);
end