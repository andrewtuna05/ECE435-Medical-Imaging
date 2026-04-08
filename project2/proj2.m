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

load("L2K.mat")

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
all_scans = reshape(BScan, N, []); 
bg_B = all_scans(:,1:D); % background columns
data_B = all_scans(:,D+1:end); % actual A-scans

Bscan_complex = generate_Ascan(bg_B, data_B, L2K);

% For both M1 and M40, they have 320 bckgrd scans so take after that
%M1_all = reshape(MScan1, N, []);
%M1_complex = generate_Ascan(M1_all(:, 1:320), M1_all(:, 321:end), L2K);
%clear M1_all MScan1;
M1_complex = single(M1_complex);

%M40_all = reshape(MScan40, N, []);
%M40_complex = generate_Ascan(M40_all(:, 1:320), M40_all(:, 321:end), L2K);
%clear M40_all MScan40;
M40_complex = single(M40_complex);

z = (0:2047) * dz; % This creates an axis from 0 to ~7mm

figure('Name', 'A-Scan Magnitude Analysis');
subplot(2,2,1);
plot(z*1e3, 20*log10(abs(Bscan_complex(:, 3000))));
title('BScan: A-Scan #3000');
xlabel('Depth (mm)'); ylabel('Magnitude (dB)'); grid on;

% --- BScan A-Scan #8000 ---
subplot(2,2,2);
plot(z*1e3, 20*log10(abs(Bscan_complex(:, 8000))));
title('BScan: A-Scan #8000');
xlabel('Depth (mm)'); ylabel('Magnitude (dB)'); grid on;

subplot(2,2,3);
plot(z*1e3, 20*log10(abs(M1_complex(:, 5000))));
title('MScan 1 A-Scan #5000');
ylabel('Magnitude (dB)'); xlabel('Depth (mm)'); grid on;

subplot(2,2,4);
plot(z*1e3, 20*log10(abs(M40_complex(:, 5000))));
title('MScan 40 A-Scan #5000');
ylabel('Magnitude (dB)'); xlabel('Depth (mm)'); grid on;

%% W/o Deconv or Background subtraction

% Extract just the one column to save memory
raw_8000 = double(data_B(:, 8000));
avg_bg = mean(double(bg_B), 2);

% Deconvolution
pixel_indices = (1:2048)';
p = polyfit(pixel_indices, avg_bg, 6);
smooth_bg = polyval(p, pixel_indices);
win = hamming(2048);

% Case 1: No BG, No Decon
% Window -> Resample -> IFFT
res_A = ifft(L2K * (raw_8000 .* win));

% Case 2: BG sub only 
% (Raw - BG) -> Window -> Resample -> IFFT
res_B = ifft(L2K * ((raw_8000 - avg_bg) .* win));

% Case 3: Full Process
res_C = ifft(L2K * (((raw_8000 - avg_bg) ./ smooth_bg) .* win));

figure('Name', 'Processing Comparison: A-Scan #8000');
hold on;

plot(z*1e3, 20*log10(abs(res_A) + eps), 'g', 'LineWidth', 1);
plot(z*1e3, 20*log10(abs(res_B) + eps), 'r', 'LineWidth', 1);                   
plot(z*1e3, 20*log10(abs(res_C) + eps), 'b', 'LineWidth', 1.5); 

title('Comparison of A-Scan Processing (A-Scan #8000)');
xlabel('Depth (mm)');
ylabel('Magnitude (dB)');
legend('No Deconv, No BG Subtraction', 'BG Subtraction Only', 'Full Process', 'Location', 'north');
grid on;
axis tight;

%% Part 3
B_processed = generate_Bscan(bg_B, data_B, L2K, true);
B_no_deconv = generate_Bscan(bg_B, data_B, L2K, false);

B_dB_proc = 20 * log10(abs(B_processed(1:1024, :)));
B_db_proc_shifted = B_dB_proc - max(B_dB_proc(:));

B_dB_none = 20 * log10(abs(B_no_deconv(1:1024, :)));
B_db_none_shifted = B_dB_none - max(B_dB_none(:));

% Gaussian filtering 
sigma = 1.2;
B_filtered = imgaussfilt(B_db_proc_shifted, sigma);

x_axis = linspace(0, bscan_width, M) * 1e3; 
y_axis = z(1:1024)* 1e3; % Depth axis matching our 1024 crop

figure('Name', 'B-Scan Comparison');

subplot(1, 2, 1);
imagesc(x_axis, y_axis, B_db_proc_shifted);
colormap gray; 
title('Unprocessed B-Scan');
xlabel('Width (mm)'); ylabel('Depth (mm)');
axis image; 

subplot(1, 2, 2);
imagesc(x_axis, y_axis, B_filtered);
colormap gray;
title('Processed B-Scan');
xlabel('Width (mm)'); ylabel('Depth (mm)');
axis image;


%% 
% For Comparison w/ and w/o Deconvolution
B_final_without = imgaussfilt(B_db_none_shifted, sigma);

figure('Name', 'Deconvolution Comparison ');

subplot(1, 2, 1);
imagesc(x_axis, y_axis, B_db_proc_shifted);
colormap gray; 
axis image;
title('B-Scan with Deconvolution');
ylabel('Depth (mm)'); xlabel('Width (mm)');
yticks(min(y_axis):0.05:max(y_axis))

subplot(1, 2, 2);
imagesc(x_axis, y_axis, B_db_none_shifted);
colormap gray; 
axis image;
title('B-Scan without Deconvolution');
ylabel('Depth (mm)'); xlabel('Width (mm)');
yticks(min(y_axis):0.2:max(y_axis))
%% Part 4

% Average A-Scan mag for M1
avg_Ascan_full = mean(abs(M1_complex), 2);
avg_Ascan = avg_Ascan_full(1:1024);
avg_Ascan_dB = 20 * log10(avg_Ascan);

figure;
plot(z(1:1024)* 1e3, avg_Ascan_dB);
title('Average A-Scan Magnitude for 1-tone M-Scan');
xlabel('Depth (mm)'); ylabel('Magnitude (dB)');
grid on;

% From figure, the highest peaks are at x = 0.216, and x = 0.3348. Since dz
% = 0.0036, this corresponds to the 61st and 94th pixel (accounting for the
% 1 index)

fs = 97656.25; % Sampling Freq

num_scans = size(M1_complex, 2);
t_raw = (0:num_scans-1) / fs;
bad_s = round(0.001 * fs); 
M_Scan_trimmed = M1_complex(:, bad_s+1:end);
t = t_raw(bad_s+1:end);

disp_1 = calc_sdpm(M_Scan_trimmed, 61, lambda_0, n);
disp_1_final = disp_1; %- disp_1(1);

L = length(disp_1_final);
Y = fft(disp_1_final);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
freq = fs*(0:(L/2))/L;
mag_fft1_dB = 20 * log10(P1);

min_freq = 0.1; % kHz
freq_khz = freq / 1000;
freq_idx_start = find(freq_khz >= min_freq, 1);

[~, rel_idx] = max(P1(freq_idx_start:end));
peak_freq1 = freq_khz(rel_idx + freq_idx_start - 1);

fprintf('Detected Speaker Tone for Pixel 61: %.2f kHz\n', peak_freq1);

figure('Name', 'MScan Analysis - Pixel 61');

subplot(3, 1, 1);
plot(t, disp_1_final, 'b-'); 
title(['Displacement at Depth 0.216 mm (Pixel 61)']);
xlabel('Time (s)'); ylabel('Displacement (nm)');
grid on;

subplot(3, 1, 2);
plot(t, disp_1_final, 'b-'); 
title('Zoomed Displacement at Depth 0.216 mm');
xlabel('Time (s)'); ylabel('Displacement (nm)');
% Zooms into the first 5ms
xlim([t(1), t(1) + 0.005]); 
grid on;

subplot(3, 1, 3);
plot(freq_khz, mag_fft1_dB, 'b-');
hold on;
plot(peak_freq1, mag_fft1_dB(rel_idx + freq_idx_start - 1), 'ro');
title('Frequency Spectrum at Depth 0.216 mm');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
xlim([min_freq, 49]); 
grid on;

disp_2 = calc_sdpm(M_Scan_trimmed, 94, lambda_0, n);
disp_2_final = disp_2 - disp_2(1);

L2 = length(disp_2_final);
Y2 = fft(disp_2_final);
P2_2 = abs(Y2/L2);
P1_2 = P2_2(1:floor(L2/2)+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);
mag_fft2_dB = 20 * log10(P1_2 + eps);

[~, rel_idx2] = max(P1_2(freq_idx_start:end));
peak_freq2 = freq_khz(rel_idx2 + freq_idx_start - 1);

fprintf('Detected Speaker Tone for Pixel 94: %.2f kHz\n', peak_freq2);

figure('Name', 'MScan Analysis - Pixel 94');

subplot(3, 1, 1);
plot(t, disp_2_final, 'r-'); 
title(['Displacement at Depth 0.334 mm (Pixel 94)']);
xlabel('Time (s)'); ylabel('Displacement (nm)');
grid on;

subplot(3, 1, 2);
plot(t, disp_2_final, 'r-'); 
title('Zoomed Displacement at Depth 0.334 mm');
xlabel('Time (s)'); ylabel('Displacement (nm)');
% Zooms into the first 5ms
xlim([t(1), t(1) + 0.005]); 
grid on;

subplot(3, 1, 3);
plot(freq_khz, mag_fft2_dB, 'r-');
hold on;
plot(peak_freq2, mag_fft2_dB(rel_idx2 + freq_idx_start - 1), 'ko');
title('Frequency Spectrum at Depth 0.334 mm');
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
xlim([min_freq, 49]); 
grid on;

%% Part 4 M40

avg_Ascan40_full = mean(abs(M40_complex), 2);
avg_Ascan40 = avg_Ascan40_full(1:1024);
avg_Ascan40_dB = 20 * log10(avg_Ascan40);

figure;
plot(z(1:1024)* 1e3, avg_Ascan40_dB);
title('Average A-Scan Magnitude for 40-tone M-Scan');
xlabel('Depth (mm)'); ylabel('Magnitude (dB)');
grid on;

bad_s = round(0.001 * fs); 
M40_trimmed = M40_complex(:, bad_s+1:end);

num_scans_40 = size(M40_complex, 2); 
t_40_raw = (0:num_scans_40-1) / fs;
t = t_40_raw(bad_s+1:end);

disp_1_40 = calc_sdpm(M40_trimmed, 61, lambda_0, n);
disp_1_40final = disp_1_40;

L_40 = length(disp_1_40final);
Y_40 = fft(disp_1_40final);
P2_40 = abs(Y_40/L_40);
P1_40 = P2_40(1:floor(L_40/2)+1);
P1_40(2:end-1) = 2*P1_40(2:end-1);
freq = fs*(0:(L_40/2))/L_40;
freq_khz = freq / 1000;
mag_fft1_dB = 20 * log10(P1_40);

min_freq = 0.1; % kHz
start_idx = find(freq_khz >= min_freq, 1);

[pks, locs] = findpeaks(P1_40(start_idx:end), 'NPeaks', 40, 'SortStr', 'descend');
unsorted_freqs = freq_khz(locs + start_idx - 1);
unsorted_mags_dB = 20 * log10(pks);

[peak1_freqs_40, sort_idx] = sort(unsorted_freqs);
peak1_mags_dB_sorted = unsorted_mags_dB(sort_idx);

figure('Name', 'MScan40 Analysis - Pixel 61');

subplot(3, 1, 1);
plot(t, disp_1_40final, 'b-'); 
title(['Displacement at Depth 0.216 mm (Pixel 61)']);
xlabel('Time (s)'); ylabel('Displacement (nm)');
grid on;

subplot(3, 1, 2);
plot(t, disp_1_40final, 'b-'); 
title('Zoomed Displacement at Depth 0.216 mm');
xlabel('Time (s)'); ylabel('Displacement (nm)');
xlim([t(1), t(1) + 0.005]); 
grid on;

subplot(3, 1, 3);
plot(freq_khz, mag_fft1_dB, 'b-');
hold on;
plot(peak1_freqs_40, peak1_mags_dB_sorted, 'ro', 'MarkerSize', 4); 
title(['Frequency Spectrum at Depth 0.216 mm (', num2str(length(peak1_freqs_40)), ' Tones)']);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
xlim([min_freq, 49]); 
grid on;

% Pixel 94 for M40 
disp_2 = calc_sdpm(conj(M40_trimmed), 94, lambda_0, n); 
disp_2_final = disp_2 - disp_2(1); 

L_40_2 = length(disp_2_final);
Y_40_2 = fft(disp_2_final);
P2_40 = abs(Y_40_2/L_40_2);
P2_40 = P2_40(1:floor(L_40_2/2)+1);
P2_40(2:end-1) = 2*P2_40(2:end-1);

freq_2 = fs * (0:(L_40_2/2)) / L_40_2;
freq_khz_2 = freq_2 / 1000;
mag_fft2_dB = 20 * log10(P2_40); 
min_freq = 0.1; 
start_idx_2 = find(freq_khz_2 >= min_freq, 1);

[pks2_dB, locs2_rel] = findpeaks(mag_fft2_dB(start_idx_2:end), 'NPeaks', 40, 'SortStr', 'descend');
locs2_final = locs2_rel + start_idx_2 - 1;
peak_freqs_40_px2 = freq_khz_2(locs2_final);

figure('Name', 'MScan Analysis - Pixel 94');

subplot(3, 1, 1);
plot(t, disp_2_final, 'r-'); 
title(['Displacement at Depth 0.3348 mm (Pixel 94)']);
xlabel('Time (s)'); ylabel('Displacement (nm)');
grid on;

subplot(3, 1, 2);
plot(t, disp_2_final, 'r-'); 
title('Zoomed Displacement at Depth 0.3348 mm');
xlabel('Time (s)'); ylabel('Displacement (nm)');
xlim([t(1), t(1) + 0.005]); 
grid on;

subplot(3, 1, 3);
plot(freq_khz_2, mag_fft2_dB, 'r-');
hold on;
plot(peak_freqs_40_px2, pks2_dB, 'ko', 'MarkerSize', 4); 

title(['Frequency Spectrum at Depth 0.3348 mm (40 Tones)']);
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
xlim([min_freq, 49]); 
grid on;



%% Functions

function A_scans_complex = generate_Ascan(bg, data, L2K)

    N = 2048;
    
    % Avg background
    t_bg = tic;
    avg_bg = mean(bg, 2); 
    data_sub = data - avg_bg;
    time_bg = toc(t_bg);

    % Deconvolution
    t_dec = tic;
    pixel_indices = (1:N)';
    poly_coeffs = polyfit(pixel_indices, avg_bg, 3);
    smooth_bg = polyval(poly_coeffs, pixel_indices);
    data_norm = data_sub ./ smooth_bg; 
    time_dec = toc(t_dec);

    % Windowing
    t_win = tic;
    win = hamming(N);
    data_win = data_norm .* win; 
    time_win = toc(t_win);
    
    % Transform to k space
    t_l2k = tic;
    data_k = L2K * data_win;
    time_l2k = toc(t_l2k);
    
    % Fourier transform from k to spatial
    t_fft = tic;
    A_scans_full = ifft(data_k, N, 1);
    A_scans_complex = A_scans_full;
    time_fft = toc(t_fft);

    % Print the timing results 
    fprintf('\n Processing Steps Time\n');
    fprintf('BG Subtraction:    %.4f seconds\n', time_bg);
    fprintf('Deconvolution:     %.4f seconds\n', time_dec);
    fprintf('Windowing:         %.4f seconds\n', time_win);
    fprintf('L2K Resampling:    %.4f seconds\n', time_l2k);
    fprintf('FFT Processing:    %.4f seconds\n', time_fft);
end

function Bscan_complex = generate_Bscan(bg, data, L2K, deconv)
    N = 2048;
    
    avg_bg = mean(bg, 2); 
    
    if deconv % Boolean tag
        % If doing deconv
        pix = (1:N)';
        p = polyfit(pix, avg_bg, 3);
        smooth_bg = polyval(p, pix);
        
        data_processed = (data - avg_bg) ./ smooth_bg;
    else
        % If not, only background subtraction
        data_processed = data - avg_bg;
    end

    data_k = L2K * (data_processed .* hamming(N));
    Bscan_complex = ifft(data_k, N, 1);
end

function d_nm = calc_sdpm(M_complex, pixel_idx, lambda0, n)
    complex_ts = M_complex(pixel_idx, :);
    phi_raw = angle(complex_ts);
    phase = unwrap(phi_raw);
    
    % Eq. 1.18 from Thesis
    d_nm = (lambda0 * phase) / (4 * pi * n); 
    
    d_nm = d_nm * 1e9; % Scale to nm
    
    d_nm = d_nm - d_nm(1); % Zero signal
end