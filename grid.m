clc; clear;

%Turn off all warnings
warning('off', 'all');

% Import the necessary libraries
import map.geodesy.*;
import matlab.net.*;
import matlab.net.http.*;

% Define the coordinates of the Fukushima Daiichi Nuclear Power Plant
fukushima_coord = [141.0333, 37.4214];

% Define the grid size in kilometers and buffer radius
grid_size_km = 100; % 约20km
buffer_km = 2000; % 1000 公里

% Calculate boundary coordinates
xmin = fukushima_coord(1) - buffer_km / 111.32;
xmax = fukushima_coord(1) + buffer_km / 111.32;
ymin = fukushima_coord(2) - buffer_km / 111.32;
ymax = fukushima_coord(2) + buffer_km / 111.32;

% Calculate the number of rows and columns of the grid
rows = round((ymax - ymin) / (grid_size_km / 111.32));
cols = round((xmax - xmin) / (grid_size_km / 111.32));
sz = [rows,cols];
% Creates the polygons of the mesh
grid_polygons = [];
for i = 1:cols
    for j = 1:rows
        x1 = xmin + (i - 1) * grid_size_km / 111.32;
        x2 = xmin + i * grid_size_km / 111.32;
        y1 = ymin + (j - 1) * grid_size_km / 111.32;
        y2 = ymin + j * grid_size_km / 111.32;
        location=[(x1+x2)/2,(y1+y2)/2];%Polygon center coordinates
        grid_polygons = [grid_polygons; struct('shape',polyshape([x1 x2 x2 x1], [y1 y1 y2 y2]),'location',location)];
    end
end

% Download and load Natural Earth data
url = 'https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip';
zip_file_path = 'ne_110m_admin_0_countries.zip';
options = weboptions('Timeout', Inf);
websave(zip_file_path, url, options);
unzip(zip_file_path, 'ne_110m_admin_0_countries');

% Read Shapefile
shapefile_path = fullfile('ne_110m_admin_0_countries', 'ne_110m_admin_0_countries.shp');
world = shaperead(shapefile_path);

% Merge the land areas of all countries into a single polygon
all_land_poly = polyshape();
for i = 1:length(world)
    current_poly = polyshape(world(i).X, world(i).Y);
    current_poly = rmholes(current_poly); % Remove holes in the polygon
    all_land_poly = union(all_land_poly, current_poly);
end

% Simplify merged polygons
all_land_poly = simplify(all_land_poly);
cj=cols*rows;
% Strips out the overlapping parts of the grid that overlap all land
sea_grid_polygons = [];
for i = 1:cj
    if ~isempty(grid_polygons(i).shape.Vertices)
        %Calculate the difference
        diff_poly = subtract(grid_polygons(i).shape, all_land_poly);
%         sea_grid_polygons = [sea_grid_polygons, diff_poly];
        sea_grid_polygons = [sea_grid_polygons; struct('diff_poly',diff_poly,'location',grid_polygons(i).location)];
        %grid_polygons = [grid_polygons; struct('shape',polyshape([x1 x2 x2 x1], [y1 y1 y2 y2]),'location',location)];
    end
end
%  Draw a map
h=figure;
hold on;
mapshow(world, 'DisplayType', 'polygon', 'FaceColor', 'white', 'EdgeColor', 'black');
for i = 1:cj
    plot(sea_grid_polygons(i).diff_poly, 'FaceColor', 'none', 'EdgeColor', [0.5, 0.5, 0.5]); % Gray grid
end
scatter(fukushima_coord(1), fukushima_coord(2), 'filled', 'r'); % The Fukushima location is marked in red
hold off;

%%Completely remove the mesh that overlaps with the land area to prepare for the simulation of cellular automata
sp_grid=[];%Marine section
ld_grid=[];%There is an overlapping grid with land
num_spgrid=0;%The number of grids in the ocean section
num_ldgrid=0;%There is a number of overlapping grids with land
for i=1:cj
    if mean(sea_grid_polygons(i).diff_poly.Vertices(:,1)) == sea_grid_polygons(i).location(1)
        %If the coordinates of the center point change after the intersection with the land, it means that the   %grid overlaps with the land and should be removed. Otherwise, it will be retained
        sp_grid=[sp_grid;sea_grid_polygons(i)];
        num_spgrid=num_spgrid+1;
    else
        ld_grid=[ld_grid;sea_grid_polygons(i)];
        num_ldgrid=num_ldgrid+1;
    end
end


save("sp_grid.mat",'sp_grid');
save("ld_grid.mat",'ld_grid');
save('num_spgrid.mat','num_spgrid');
save('num_ldgrid.mat','num_ldgrid');
save('grid_polygons.mat','grid_polygons');
saveas(h,'map','fig');
save('grid_size_km.mat','grid_size_km');
save('sz.mat','sz');
