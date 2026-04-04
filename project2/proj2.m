% Andrew Yuan ECE435 - Project 2

close all; 
clear; 
clc;

%% Data Loading

fileID_Bscan = fopen("BScan_Layers.raw");
fileID_Mscan1 = fopen("MScan1.raw");
fileID_Mscan40 = fopen("MScan40.raw");

BScan = fread(fileID_Bscan, 'uint16');
MScan1 = fread(fileID_Mscan1, 'uint16');
MScan40 = fread(fileID_Mscan40, 'uint16');


%% Part 1

% Parameters
dz = 3.6e-6; % Axial pixel size 
NA = 0.055; % Numerical Aperture
lambda_0 = 1300e-9; % Center Wavelength
d_lambda = 100e-9; % Bandwidth
n = 1; % Refractive Index of Air
C = 0.37;

% Lateral resolution
lat_res = C*lambda_0/NA; % C = 0.37
% The lateral resolution is limited by the numerical aperture and center
% wavelength of the imaging system. We want higher resolution (more detailed) so we would want higher
% numerical aperture (which results in finer more detailed images) and a
% smaller center wavelength. The tradeoff is that this lowers our penetration depth.

% Axial resolution
ax_res = 2*log(2)*(lambda_0).^2/(n*pi*d_lambda);
% The axial resolution is limited by the source bandwidth. For a given
% center wavelength, having a larger bandwidth would decrease the axial
% resolution, corresponding to higher (better) resolution. 

% B-Scan pixel aspect ratios
D = 175; % # of backgrounds
N = 2048; % # of pixels per line
bscan_width = 0.001;
[total , ~] = size(BScan);

M = total/2048 - D; % # of A-scans

total_axial_dep = dz * 2048;
total_axial_dep_mir = total_axial_dep/2;
lateral_pixel_size = (bscan_width/M);

pixel_aspect_ratio = dz/lateral_pixel_size;

fprintf("Pixel Aspect Ratio: %.2f\n", pixel_aspect_ratio);
fprintf("B-Scan Aspect Ratio: %.3f x %.6f\n", bscan_width, total_axial_dep);
fprintf("B-Scan Aspect (Mirror) Ratio: %.3f x %.6f\n", bscan_width, total_axial_dep_mir);

%% Part 2

% Manipulating Shape of Ascan



function gen_complex_Ascan = 