clc
close all
clear all

Img = imread("Image1.tiff");
res = 30; %LANDSAT-8 SWIR-2 30m

GS_Ratio = 4; % Ratio of Ground Swath of HRR/HAR
X_cut = 35; % Ground Swath Required [km]
Y_cut = 100; % Azimuth Length of Image Required [km]

% The HiRes of both Antenna = res of Image (30m)
HRR_AR = 4; % Ratio of LowRes/HiRes in HRR
HAR_AR = 5; % Ratio of LowRes/HiRes in HAR

%% Cropping Image
Img_size = size(Img);
Img_size_km = (Img_size.*res)./1e3;

Img_cropped_size = ceil([X_cut*1e3/res Y_cut*1e3/res]); % Cropped Image Size in pixels
Img_cent = ceil(Img_size./2); % Image Center in Pixels
crop_TopLeft = ceil(Img_cent - (Img_cropped_size./2)); % TopLeft Pixel location of Cropping Rectangle
cropping_rect = [crop_TopLeft Img_cropped_size]; % Cropping Rectangle Struct

Img_cropped = imcrop(Img, cropping_rect);
Img_cropped_size = size(Img_cropped);

rems = [rem(Img_cropped_size(1),HRR_AR) rem(Img_cropped_size(2),HAR_AR)]; % Cropping extra pixels to make
% it cropped image compatible to the HRR_AR and HAR_AR values.

Img_cropped = Img_cropped(1:(end-rems(1)), 1:(end-rems(2)));
Img_cropped_size = size(Img_cropped);


%% Image normalisation
im = double(Img_cropped);
Img_cropped = (im-min(im(:)))/(max(im(:))-min(im(:)));

%% Simulate Acquisition (Pixel Average)
% Simulate HAR Image
% Combining Columns
% No. of Columns in cropped image Must be Divisible by HAR_AR

Img_HAR_size = [Img_cropped_size(1) (Img_cropped_size(2)/HAR_AR)]; % HRR Image size in Pixels
Img_HAR_mean = zeros(Img_HAR_size);
for i=1:Img_HAR_size(2)
    k = HAR_AR*(i-1);
    a = k+1;
    b = k+HAR_AR;
    Img_HAR_mean(:,i) = mean(Img_cropped(:,a:b),2); % Pixel Average Method
    %Img_HAR_sum(:,i) = uint16(sum(Img_cropped(:,a:b),2)); % Pixel Sum Method
end

% Simulate HRR Image
% Combining Rows
% No. of Rows in cropped image Must be Divisible by HAR_AR

Img_HRR_size = [(Img_cropped_size(1)/HRR_AR) Img_cropped_size(2)]; % HRR Image size in Pixels
Img_HRR_mean = zeros(Img_HRR_size);
for i=1:Img_HRR_size(1)
    k = HRR_AR*(i-1);
    a = k+1;
    b = k+HRR_AR;
    Img_HRR_mean(i,:) = mean(Img_cropped(a:b,:),1); % Pixel Average Method
    %Img_HRR_sum(i,:) = uint16(sum(Img_cropped(a:b,:),1)); % Pixel Sum Method
end

%% Simulate SAR DEM Generation
% Dummy HAR Image
Dummy_HAR_Img = zeros(Img_cropped_size);
for i = 1:Img_HAR_size(2)
    k = HAR_AR*(i-1);
    a = k+1;
    b = k+HAR_AR;
    for j=a:b
        Dummy_HAR_Img(:,j) = Img_HAR_mean(:,i);
    end
end

% Dummy HRR Image
Dummy_HRR_Img = zeros(Img_cropped_size);
for i = 1:Img_HRR_size(1)
    k = HRR_AR*(i-1);
    a = k+1;
    b = k+HRR_AR;
    for j=a:b
        Dummy_HRR_Img(j,:) = Img_HRR_mean(i,:);
    end
end

% Mean Dummy
Constructed_Img = (Dummy_HAR_Img + Dummy_HRR_Img)./2;

%% Error
error = (Img_cropped - Constructed_Img);
error_mean = mean(double(error),'all');
error_var = var(double(error),0,'all');
error_std = std(double(error),0,'all');
disp('Error b/w Ground Truth & Generated Images:')
disp('    Mean     Variance   Standard Deviation')
disp([error_mean error_var error_std])

%% Single scan Simulated Image
cropped_cent = ceil(Img_cropped_size./2); % Center of Cropped Image

HAR_cent = ceil(Img_HAR_size/2);

SS_HAR_width = ceil(Img_HAR_size(2)/GS_Ratio); % Width of HAR Scan Ground Swath in Pixels
SS_HAR_length = Img_HAR_size(1);
SS_HAR_size = [SS_HAR_length SS_HAR_width];

SS_HAR_crop_rect = ceil([fliplr(HAR_cent-(SS_HAR_size./2)) fliplr(SS_HAR_size)]);

Img_SS_HAR = imcrop(Img_HAR_mean,SS_HAR_crop_rect);

Img_SS_HRR = Img_HRR_mean;

%% Plot Results
figure(1)

subplot(1,2,1),
imshow(Img)
hold on;
rectangle('Position', cropping_rect, 'EdgeColor','y');
hold off;
title('Original Image')

subplot(1,2,2),
imshow(Img_cropped)
title('Cropped Image')

T1 = {'Reference Images - (' num2str(res) 'm)'};
sgtitle(T1)

figure(2)

subplot(1,3,1),
imshow(Img_HAR_mean)
title('HAR Scan - Complete Ground Swath')
hold on;
rectangle('Position', SS_HAR_crop_rect, 'EdgeColor','red');
hold off;

subplot(1,3,2),
imshow(Img_SS_HAR)
title('Single Scan HAR Image')

subplot(1,3,3),
imshow(Img_SS_HRR)
title('Single Scan HRR Image')

T1 = {'Simulated Single Scan Images'};
sgtitle(T1)

figure(3)

subplot(1,3,1),
imagesc(Img_HAR_mean)
title('HAR Scan - Complete Ground Swath')
hold on;
rectangle('Position', SS_HAR_crop_rect, 'EdgeColor','red');
hold off;

subplot(1,3,2),
imagesc(Img_SS_HAR)
title('Single Scan HAR Image')

subplot(1,3,3),
imagesc(Img_SS_HRR)
title('Single Scan HRR Image')

T1 = {'Simulated Single Scan Images - Equal Pixel Size'};
sgtitle(T1)

figure(4)
subplot(1,2,1),
imshow(Img_HAR_mean)
title('HAR Image')

subplot(1,2,2),
imshow(Img_HRR_mean)
title('HRR Image')

T1 = {'Simulated Acquisitions of Complete Ground Swath'};
sgtitle(T1)

figure(5)
subplot(1,3,1),
imshow(Dummy_HAR_Img)
title('Dummy HAR Image')

subplot(1,3,2),
imshow(Dummy_HRR_Img)
title('Dummy HRR Image')

subplot(1,3,3),
imshow(Constructed_Img)
title('Constructed Image')

T1 = {'Constructed Synthetic Image of Complete Ground Swath'};
sgtitle(T1)

figure(6)
subplot(1,3,1),
imshow(Img_cropped)
title('Cropped Image')

subplot(1,3,2),
imshow(Constructed_Img)
title('Constructed Image')

subplot(1,3,3),
imshow(error)
title('Absolute Error')
T1 = {'Constructed Synthetic Image Comparission'};
sgtitle(T1)

