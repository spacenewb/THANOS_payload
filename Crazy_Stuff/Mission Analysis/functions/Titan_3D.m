function Titan_3D(image_in)
% Plot Titan 3-D
%
% Inputs:
% image_in      : image of Titan to plot over the sphere
%
% Plots:
% Titan 3-D
% 
% Contributors:
% Gaballo Paolo
%  
% Versions:
% 2021-04-12, second version
% 2021-02-13, first version (earth)


if nargin == 0
    image_in = 'TitanTexture.jpg';
end

Rt = 2574.73;  

% Define the number of panels to be used to model the sphere 
npanels = 180;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, npanels);

% Create the globe with the surf function
globe = surf(x, y, -z, 'EdgeColor', 'none', 'FaceColor', 'none');
hold on;
axis equal;
%Earth_image = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

cdata = imread(image_in);
% cdata = imread('EarthTexture.jpg');

set(globe,'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha',1);
end