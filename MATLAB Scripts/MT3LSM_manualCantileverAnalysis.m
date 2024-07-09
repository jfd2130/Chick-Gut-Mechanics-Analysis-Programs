function [calcData,rawData] = MT3LSM_manualCantileverAnalysis
%% Setup

% I. allow the user to select an image file (start of run, usually image 3):
[fileName, pathName] = uigetfile('*.jpg', 'Open first image of run:');
fullName = fullfile(pathName,fileName); % these four lines are from Meshoo on MATLAB Answers...
filelist = dir([fileparts(fullName) filesep '*.jpg']);
fileNames = {filelist.name}';
fileNames = fileNames(1:1:length(fileNames));
images = numel(filelist);
warning('off','images:initSize:adjustingMag'); % reduces Command Window clutter

% VI. initialize a 3D array (m) to hold measured values:
rawData = zeros(images,3,2); % dimensions (in order): image number, measurement group (top, bottom, tissue), point within group (L,R)


%% Initialization:

resolution = 0.05; % mm/px

% I. cantilever parameters for the run:
R = 0.064*10^-3; % get the radius (R) as a number (m)
L = 30*10^-3; % get the length (L) as a number (m)
E = 410e9; % known modulus (E) of tungsten (Pa)
% use R, L, and E to calculate bending stiffness:
k_b = (3*E*pi*R^4)/(4*L^3);


%% Image Processing

% select points on cantilevers and tissue:
topPoints = zeros(images,2,2); % image, point (L,R), coordinate (x,y)
bottomPoints = zeros(images,2,2); % image, point (L,R), coordinate (x,y)
tissuePoints = zeros(images,2,2); % image, point (L,R), coordinate (x,y)
figure(1)
imgBf = imread(fullfile(pathName,fileNames{1}));
imgG = imcomplement(imgBf(:,:,2));
imshow(imgG); hold on
for i = 1:2 % manually select top points (L then R)
    topPoints(1,i,:) = uint16(ginput(1));
    scatter(topPoints(1,i,1),topPoints(1,i,2),60,'r','+');
end
for i = 1:2 % manually select bottom points (L then R)
    bottomPoints(1,i,:) = uint16(ginput(1));
    scatter(bottomPoints(1,i,1),bottomPoints(1,i,2),60,'b','+');
end
for i = 1:2 % manually select tissue points (L then R)
    tissuePoints(1,i,:) = uint16(ginput(1));
    scatter(tissuePoints(1,i,1),tissuePoints(1,i,2),25,'g','o');
end
pause(0.5);
delete(gca);
% save ROIs around selected points:
dim_ROI = 10;
topPointL_ROI = imgG(topPoints(1,1,2)-dim_ROI:topPoints(1,1,2)+dim_ROI, topPoints(1,1,1)-dim_ROI:topPoints(1,1,1)+dim_ROI);
topPointR_ROI = imgG(topPoints(1,2,2)-dim_ROI:topPoints(1,2,2)+dim_ROI, topPoints(1,2,1)-dim_ROI:topPoints(1,2,1)+dim_ROI);
bottomPointL_ROI = imgG(bottomPoints(1,1,2)-dim_ROI:bottomPoints(1,1,2)+dim_ROI, bottomPoints(1,1,1)-dim_ROI:bottomPoints(1,1,1)+dim_ROI);
bottomPointR_ROI = imgG(bottomPoints(1,2,2)-dim_ROI:bottomPoints(1,2,2)+dim_ROI, bottomPoints(1,2,1)-dim_ROI:bottomPoints(1,2,1)+dim_ROI);
tissuePointL_ROI = imgG(tissuePoints(1,1,2)-dim_ROI:tissuePoints(1,1,2)+dim_ROI, tissuePoints(1,1,1)-dim_ROI:tissuePoints(1,1,1)+dim_ROI);
tissuePointR_ROI = imgG(tissuePoints(1,2,2)-dim_ROI:tissuePoints(1,2,2)+dim_ROI, tissuePoints(1,2,1)-dim_ROI:tissuePoints(1,2,1)+dim_ROI);
% subplot(3,2,1); imshow(histeq(topPointL_ROI))
% subplot(3,2,2); imshow(histeq(topPointR_ROI))
% subplot(3,2,3); imshow(histeq(bottomPointL_ROI))
% subplot(3,2,4); imshow(histeq(bottomPointR_ROI))
% subplot(3,2,5); imshow(histeq(tissuePointL_ROI))
% subplot(3,2,6); imshow(histeq(tissuePointR_ROI))

for i = 2:images
    % 1. import image:
    imgBf = imread(fullfile(pathName,fileNames{i}));
    
    % 2. split RGB channels (just keep G):
    imgG = imcomplement(imgBf(:,:,2));

    imshow(imgG); hold on
    scatter(topPoints(i-1,:,1),topPoints(i-1,:,2),30,'r','+');
    scatter(bottomPoints(i-1,:,1),bottomPoints(i-1,:,2),30,'b','+');
    scatter(tissuePoints(i-1,:,1),tissuePoints(i-1,:,2),12,'g','o');
    for j = 1:2 % manually select top points (L then R)
        topPoints(i,j,:) = uint16(ginput(1));
        scatter(topPoints(i,j,1),topPoints(i,j,2),60,'r','+');
    end
    for j = 1:2 % manually select bottom points (L then R)
        bottomPoints(i,j,:) = uint16(ginput(1));
        scatter(bottomPoints(i,j,1),bottomPoints(i,j,2),60,'b','+');
    end
    for j = 1:2 % manually select tissue points (L then R)
        tissuePoints(i,j,:) = uint16(ginput(1));
        scatter(tissuePoints(i,j,1),tissuePoints(i,j,2),25,'g','o');
    end
    pause(0.01);
    delete(gca);
end

% 4. save x values as distance measurments (mm units):
rawData(:,1,:) = topPoints(:,:,1);
rawData(:,2,:) = bottomPoints(:,:,1);
rawData(:,3,:) = tissuePoints(:,:,1);
rawData = rawData .* resolution;

figure(1)
imgBf = imread(fullfile(pathName,fileNames{1}));
imshow(imcomplement(imgBf(:, :, 2))); hold on
scatter(topPoints(:,:,1),topPoints(:,:,2),30,'r','+');
scatter(bottomPoints(:,:,1),bottomPoints(:,:,2),30,'b','+');
scatter(tissuePoints(:,:,1),tissuePoints(:,:,2),13,'g','o');

figure(2)
imgBf = imread(fullfile(pathName,fileNames{images}));
imshow(imcomplement(imgBf(:, :, 2))); hold on
scatter(topPoints(:,:,1),topPoints(:,:,2),30,'r','+');
scatter(bottomPoints(:,:,1),bottomPoints(:,:,2),30,'b','+');
scatter(tissuePoints(:,:,1),tissuePoints(:,:,2),13,'g','o');

figure(3)
timePoints = linspace(0,1,length(topPoints(:,:,1)));
plot(timePoints,topPoints(:,1,1)-topPoints(1,1,1),':r','Linewidth',6); hold on
plot(timePoints,topPoints(:,2,1)-topPoints(1,2,1),'r','Linewidth',4);
plot(timePoints,bottomPoints(:,1,1)-bottomPoints(1,1,1),':b','Linewidth',4);
plot(timePoints,bottomPoints(:,2,1)-bottomPoints(1,2,1),'b','Linewidth',2);
plot(timePoints,tissuePoints(:,1,1)-tissuePoints(1,1,1),':g','Linewidth',4);
plot(timePoints,tissuePoints(:,2,1)-tissuePoints(1,2,1),'g','Linewidth',2);


%% Calculations:

% stretch and stress:
distance_L_initial = abs(rawData(1,1,1)-rawData(1,2,1));
distance_R_initial = abs(rawData(1,1,2)-rawData(1,2,2));
displacement_L = abs(rawData(:,1,1)-rawData(:,2,1)) - distance_L_initial;
displacement_R = abs(rawData(:,1,2)-rawData(:,2,2)) - distance_R_initial;
displacement = (displacement_L + displacement_R) ./ 2; % use the mean of the two cantilevers
length_initial = abs(rawData(1,3,1)-rawData(1,3,2));
length_current = abs(rawData(:,3,1)-rawData(:,3,2));
% Set increasing displacement constraint, then set new initial displacement to 0:
% valRange = 5;
% for i_image = 1:length(displacement)
%     fit = fitlm(linspace(0,1,valRange)', displacement(i_image:i_image+valRange-1), 'poly1'); % linear regression
%     if fit.Coefficients.Estimate(2) > 0 && fit.Coefficients.pValue(2) < 0.05 % i.e. the slope is confidently (95%) positive in a window of "valRange" number of measurements
%         % set measurement i_image as the first measurement (when displacement = 0):
%         distance_L_initial = abs(rawData(i_image,1,1)-rawData(i_image,2,1));
%         distance_R_initial = abs(rawData(i_image,1,2)-rawData(i_image,2,2));
%         displacement_L = abs(rawData(i_image:end,1,1)-rawData(i_image:end,2,1)) - distance_L_initial;
%         displacement_R = abs(rawData(i_image:end,1,2)-rawData(i_image:end,2,2)) - distance_R_initial;
%         displacement = (displacement_L + displacement_R) ./ 2; % use the mean of the two cantilevers
%         length_initial = abs(rawData(i_image,3,1)-rawData(i_image,3,2));
%         length_current = abs(rawData(i_image:end,3,1)-rawData(i_image:end,3,2));
% 
%         disp([ 'starting image: ' num2str(i_image) ]); % output the starting image
%         break
%     end
% end
% compute final values:
stretch = length_current ./ length_initial;
force = displacement .* k_b;
area = 0.12 * 0.82; % estimated as: mean DMSO sample thickness * mean DMSO sample width
stress = force ./ area;

% output:
calcData = [stretch stress];


end