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

% squeeze to make 2-D image from slice
% img = squeeze(stack(:,:, 1));
% imshow(img)


%% Part 2
sample_img = im2gray(squeeze(stack(:,:, 124)));

% 3x3 Median Filter
med_fil_3 = medfilt2(sample_img);
figure;
imshow(med_fil_3);
title('Median Filtered 3x3');

% 5x5 Median Filter
med_fil_5 = medfilt2(sample_img, [5 5]);
figure;
imshow(med_fil_5);
title('Median Filtered 5x5');

% 3x3 Mean Filter
mean_fil_3 = imboxfilt(sample_img);
figure;
imshow(mean_fil_3);
title('Mean Filtered 3x3');

% 5x5 Mean Filter 
mean_fil_5 = imboxfilt(sample_img, 5);
figure;
imshow(mean_fil_5);
title('Mean Filtered 5x5');

% Salt-n-Pepper noise
sp_img = zeros( 512, 512, 3,'like', sample_img); % Creating 3 image templates 

p_vals = [0.01, 0.05, 0.15];

for i = 1:length(p_vals)
    p = p_vals(i);
    sp_img(:, :, i) = imnoise(sample_img, "salt & pepper", 2*p); % I think imnoise splits between white and black evenly so it works here
    figure;
    imshow(sp_img(:, :, i))
    title(sprintf('Salt-n-Pepper Noise p = %.2f', p))
end

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
