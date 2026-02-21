% Andrew Yuan ECE435 - Project 1

close all; 
clear; 
clc;

% Loading images
imds = imageDatastore("thoraxCT", 'FileExtensions', '.jpg'); % nice it auto-sorted!

%% Part 1
img = im2gray(readimage(imds, 1));
stack = zeros(512, 512, 237, 'like', img);

for i = 1:237 
    new_img = im2gray(readimage(imds, i));
    stack(:, :, i) = new_img;
end
disp(size(stack)); % 512 x 512 x 237 Good!

dx = 0.703;
dy = 0.703;
dz = 0.625;

img1 = squeeze(stack(215,:,:));
figure;
imshow(img1)
axis on;

%squeeze to make 2-D image from slice
figure;
subplot(1,3,1);
img1 = squeeze(stack(215,:,:));
imshow(img1, 'XData', [0 236*dx],'YData', [0 511*dy])
axis on
axis image
xlabel('z (mm)')
ylabel('y (mm)')
title('x = 215');

subplot(1,3,2);
img2 = squeeze(stack(:,215,:));
imshow(img2, 'XData', [0 236*dx],'YData', [0 511*dy])
axis on
axis image
xlabel('z (mm)')
ylabel('x (mm)')
title('y = 215');

subplot(1,3,3);
img3 = squeeze(stack(:,:,130));
imshow(img3, 'XData', [0 511*dx],'YData', [0 511*dy])
axis on
axis image
xlabel('x (mm)')
ylabel('y (mm)')
title('z = 130');

%% Part 2
sample_img = im2gray(squeeze(stack(:,:, 130)));

figure;

% 3x3 Median Filter
subplot(1,2,1);
med_fil_3 = medfilt2(sample_img);
imshow(med_fil_3, []);
title('Median Filtered 3x3');

% 3x3 Mean Filter
subplot(1,2,2);
mean_fil_3 = imboxfilt(sample_img);
imshow(mean_fil_3, []);
title('Mean Filtered 3x3');

figure;

% 5x5 Median Filter
subplot(1,2,1);
med_fil_5 = medfilt2(sample_img, [5 5]);
imshow(med_fil_5, []);
title('Median Filtered 5x5');

subplot(1,2,2);
mean_fil_5 = imboxfilt(sample_img, 5);
imshow(mean_fil_5, []);
title('Mean Filtered 5x5');


% Salt-n-Pepper noise
sp_img = zeros( 512, 512, 3,'like', sample_img); % Creating 3 image templates 

p_vals = [0.01, 0.05, 0.15];

for i = 1:length(p_vals)
    p = p_vals(i);
    sp_img(:, :, i) = imnoise(sample_img, "salt & pepper", 2*p); % I think imnoise splits between white and black evenly so it works here
end

figure;
imshow(sp_img(:, :, 1))
title(sprintf('Salt-n-Pepper Noise p = 0.01'));

figure;
imshow(sp_img(:, :, 2));
title(sprintf('Salt-n-Pepper Noise p = 0.05'));

figure;
imshow(sp_img(:, :, 3));
title(sprintf('Salt-n-Pepper Noise p = 0.15'));


% Applying filters to Salt-n-Pepper images
sp1_med_fil_3 = medfilt2(sp_img(:, :, 1));

sp2_med_fil_3 = medfilt2(sp_img(:, :, 2));

sp3_med_fil_3 = medfilt2(sp_img(:, :, 3));

sp1_mean_fil_3 = imboxfilt(sp_img(:, :, 1));

sp2_mean_fil_3 = imboxfilt(sp_img(:, :, 2));

sp3_mean_fil_3 = imboxfilt(sp_img(:, :, 3));

figure;
imshow(sp1_med_fil_3);
title('Median Filtered SnP Noise p = 0.01');

figure;
imshow(sp2_med_fil_3)
title('Median Filtered SnP Noise p = 0.05');

figure;
imshow(sp3_med_fil_3)
title('Median Filtered SnP Noise p = 0.15');

figure;
imshow(sp1_mean_fil_3)
title('Mean Filtered SnP Noise p = 0.01');

figure;
imshow(sp2_mean_fil_3)
title('Mean Filtered SnP Noise p = 0.05');

figure;
imshow(sp3_mean_fil_3)
title('Mean Filtered SnP Noise p = 0.15');
%% Part 3

gauss_img = zeros(512, 512, 2,'like', sample_img);
gauss_sample_img = double(sample_img);
% how does the gaussian noise func apply noise? like is it adding to the
% intensity value?

gauss_var = [5, 15];

for i = 1:length(gauss_var)
    k = gauss_var(i);
    gauss_noise = normrnd(0, k, [512, 512]);
    gauss_img(:, :, i) = gauss_sample_img + gauss_noise;
    if i == 1
        figure;
    end
    subplot(1,2,i);
    imshow(gauss_img(:, :, i))
    title('Gaussian Noise \sigma = %.2f', k)
end

% sigma = 0.5 --> 3 pixel gaussian filter
filt_gauss_img = imgaussfilt(gauss_img(:, :, 1), 0.5);
figure;
subplot(1,2,1)
imshow(filt_gauss_img)
title('Gaussian Noise = 5 + Gaussian Filter \sigma = 0.5 pixels');

filt_gauss_img_2 = imgaussfilt(gauss_img(:, :, 2), 1);
subplot(1,2,2)
imshow(filt_gauss_img_2)
title('Gaussian Noise = 15 + Gaussian Filter \sigma = 1 pixels');

delta_x = 2;
delta_z = 1;

sigma_xy = delta_x / (2*sqrt(2*log(2)));
sigma_z = delta_z / (2*sqrt(2*log(2)));
sigma_x_px = sigma_xy / dx; 
sigma_y_px = sigma_xy / dy;
sigma_z_px = sigma_z / dz;

img1_blur = imgaussfilt(img1, [sigma_y_px, sigma_z_px]);
img2_blur = imgaussfilt(img2, [sigma_x_px, sigma_z_px]);
img3_blur = imgaussfilt(img3, [sigma_y_px, sigma_x_px]);

figure;
subplot(1,2,1)
imshow(img1);
title(sprintf('x = 215 Original'))
subplot(1,2,2)
imshow(img1_blur);
title(sprintf('x = 215 Gaussian Blur'))

figure;
subplot(1,2,1)
imshow(img2);
title(sprintf('y = 215 Original'))
subplot(1,2,2)
imshow(img2_blur);
title(sprintf('y = 215 Gaussian Blur'))

figure;
subplot(1,2,1)
imshow(img3);
title(sprintf('z = 130 Original'))
subplot(1,2,2)
imshow(img3_blur);
title(sprintf('z = 130 Gaussian Blur'))
%% Part 4
% Regular Canny Edge
foramina_img = im2gray(squeeze(stack(:,180,:)));
figure;
subplot(1,2,1)
imshow(foramina_img);
subplot(1,2,2);
edge_foramina_img = edge(foramina_img,'Canny');
imshow(edge_foramina_img);
sgtitle('y = 180, Canny Edge Detection')

% The issue is that its not very clear lol. The outline of the foramina is
% not visible, I can only make out the foramina via the outline of the
% vertebrae and looking in between.

% 3x3 Median Filter
median_filt_foramina = medfilt2(foramina_img);
figure;
subplot(1,2,1);
imshow(median_filt_foramina);
subplot(1,2,2);
edge_median_foramina = edge(median_filt_foramina, 'Canny');
imshow(edge_median_foramina);
sgtitle('3x3 Median Filter Canny Edge Detection');

% 2-D Gaussian Filter
foramina_2D = imgaussfilt(foramina_img, [sigma_x_px, sigma_z_px]);
figure;
subplot(1,2,1);
imshow(foramina_2D);
subplot(1,2,2);
edge_foramina_2D = edge(foramina_2D, 'Canny');
imshow(edge_foramina_2D);
sgtitle('Fancy 2D Gaussian Filter Canny Edge Detection');


% Gaussian Noise
foramina_img_double = double(foramina_img);
gauss_noise = normrnd(0, 15, [512, 237]);
gauss_noise_foramina = foramina_img_double + gauss_noise;
gauss_noise_foramina = max(gauss_noise_foramina, 0);
gauss_noise_foramina = min(gauss_noise_foramina, 255);
gauss_noise_foramina = uint8(gauss_noise_foramina);
figure;
subplot(1,2,1);
imshow(gauss_noise_foramina);
subplot(1,2,2);
edge_median_foramina = edge(gauss_noise_foramina, 'Canny');
imshow(edge_median_foramina);
sgtitle('Gaussian Noise \sigma = 15 Canny Edge Detection');


% Gaussian Filter 3 pixel std.dev
figure;
subplot(1,2,1);
foramina_blur = imgaussfilt(foramina_img, [3, 3]);
imshow(foramina_blur);
subplot(1,2,2);
edge_foramina_blur = edge(foramina_blur, 'Canny');
imshow(edge_foramina_blur)
sgtitle('Gaussian Filter \sigma = 3 Canny Edge Detection');

% Regular Canny Edge to Z slice
foramina_img_z = im2gray(squeeze(stack(:,:,130)));
figure;
subplot(1,2,1)
imshow(foramina_img_z);
subplot(1,2,2);
edge_foramina_img_z = edge(foramina_img_z,'Canny');
imshow(edge_foramina_img_z);
sgtitle('z = 130, Canny Edge Detection')