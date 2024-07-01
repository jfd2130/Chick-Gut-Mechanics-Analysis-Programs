function [calcData,rawData,runInfo] = T4Axio_mesenteryNoMarkers_TEST()
%% Setup

% I. allow the user to select an image file (start of run, usually image 3):
[fileName, pathName] = uigetfile('*.tif', 'Open first image of run:');
fullName = fullfile(pathName,fileName); % these four lines are from Meshoo on MATLAB Answers...
filelist = dir([fileparts(fullName) filesep '*.tif']);
fileNames = {filelist.name}';
fileNames = fileNames(1:1:length(fileNames));
imagesAll = length(fileNames);
warning('off','images:initSize:adjustingMag'); % reduces Command Window clutter

% II. use the file name (name format: "run#_A-(sample1)_B-(sample2)_C-(sample3)_tXX.tif")
% of the selected file to find sample names and numbers:
close all
underscores = strfind(fileName,'_'); % underscores separate all elements in the filename (may also have condition!)
sA = fileName(underscores(1)+3:underscores(2)-1); % uses the above indices to find the name of the sample on cantilever A
sB = fileName(underscores(2)+3:underscores(3)-1);
sC = fileName(underscores(3)+3:underscores(4)-1);

% III. find the image number of the selected file to use for image indexing:
firstImg = 1; % stored for future reference
images = imagesAll - firstImg + 1; % to determine dimension length

% IV. set image analysis parameters:
length_ROIradius = 25;
length_ROItranslation = 15;

% V. set stringency of positive slope test for finding start of run:
valRange = 5; % X number of measurements in a row to be tested for (+) slope

% VI. initialize a 3D array (M) to hold measured values:
M = zeros(3,4,images,2); % dimensions (in order): sample (A,B,C), measurement (actuator, cantilever, tissue 1, tissue 2), image number, coordinate (x, y)
img_runStart = zeros(3,1); % for storing the image determined as the start of the run for each samples

% VII. welcome message:
disp('Welcome!');

checkImageEvery = 10;
brightenImgBy = 50;



%% Initialization:

% I. present dialogue box to enter cantilever and velocity parameters for the run:
R = 0.064*10^-3; % get the radius (R) as a number (m)
L = 41.2*10^-3; % get the length (L) as a number (m)
E = 410e9; % known modulus (E) of tungsten (Pa)
% use R, L, and E to calculate bending stiffness:
k_b = (3*E*pi*R^4)/(4*L^3);

% II. get the resolution:
res = 0.091; % resolution (mm/px)

% III. present dialogue box to enter sample radii:
thickness = 0.12; % get the thickness (t) of samples on A as a number (mm)

% IV. update message:
disp('Initialization complete. Image processing underway...');



%% Image Processing

% take user inputs to select data points on the first image:
imgBf = imread(fullfile(pathName,fileNames{1}));
imshow(imgBf); hold on
% select the actuator point:
title('Select the base:');
M(1,1,1,:) = round(uint16(ginput(1)));
M(2,1,1,:) = M(1,1,1,:);
M(3,1,1,:) = M(1,1,1,:);
scatter(M(1,1,1,1),M(1,1,1,2),60,'g','+');
% select cantilever points:
title('Select the cantilevers (left to right):');
for i = 1:3
    M(i,2,1,:) = round(uint16(ginput(1)));
    scatter(M(i,2,1,1),M(i,2,1,2),60,'g','+');
end
% select tissue points:
title('Select the top and bottom of each sample (left to right):');
for i = 1:2
    M(1,2+i,1,:) = round(uint16(ginput(1)));
    scatter(M(1,2+i,1,1),M(1,2+i,1,2),60,'w','+');
end
for i = 1:2
    M(2,2+i,1,:) = round(uint16(ginput(1)));
    scatter(M(2,2+i,1,1),M(2,2+i,1,2),60,'y','+');
end
for i = 1:2
    M(3,2+i,1,:) = round(uint16(ginput(1)));
    scatter(M(3,2+i,1,1),M(3,2+i,1,2),60,'c','+');
end
pause(0.2);
delete(gca);

% iterate through rest of images and analyse automatically:
figure(1) % to display point measurement example images as processing progresses
for i_img = 2:images
    % import images in 3 channels 
    imgBf = imread(fullfile(pathName,fileNames{i_img}));
    
    % iterate through samples:
    for i_sample = 1:size(M,1)
        % iterate through data points for the current sample:
        for i_point = 1:size(M,2)
            % save measured values:
            M(i_sample,i_point,i_img,:) = ShiftROI_ssim(i_sample,i_point,i_img,length_ROIradius);
        end
    end
    
    % for visual confirmation, show example images with measured points:
    if rem(i_img,checkImageEvery) == 0
        if (i_img)/checkImageEvery+1 <= 5*6
        % dsplay image:
        subplot(5,6,(i_img)/checkImageEvery);
        imshow(imgBf+brightenImgBy); hold on
        % show actuator:
        scatter(M(1,1,i_img,1),M(1,1,i_img,2),50,'g','+');
        % show cantilevers:
        scatter(M(1,2,i_img,1),M(1,2,i_img,2),50,'g','+');
        scatter(M(2,2,i_img,1),M(2,2,i_img,2),50,'g','+');
        scatter(M(3,2,i_img,1),M(3,2,i_img,2),50,'g','+');
        % show tissue points:
        scatter(M(1,3:4,i_img,1),M(1,3:4,i_img,2),50,'w','+');
        scatter(M(2,3:4,i_img,1),M(2,3:4,i_img,2),50,'y','+');
        scatter(M(3,3:4,i_img,1),M(3,3:4,i_img,2),50,'c','+');
        end

        disp(['image complete: ' num2str(i_img)]);
    end

end

% discount data from images before motion is detectable on actuator:
img_firstMotion = 1;
for i_img = 2:round(images/2)
    if M(1,1,i_img,2) >= M(1,1,i_img-1,2)
        img_firstMotion = i_img;
    end
end

% convert y values to distance measurments (mm units):
rawData = squeeze(M(:,:,:,2)) * res;

% update message:
disp('Measurements complete. Calculations underway...');



%% Calculations:

% Cantilever A:
[sA_lam, sA_P, img_runStart_A] = CalcStretchStress(rawData, 1);
img_runStart(1) = img_runStart_A;

% Cantilever B:
[sB_lam, sB_P, img_runStart_B] = CalcStretchStress(rawData, 2);
img_runStart(2) = img_runStart_B;

% Cantilever C:
[sC_lam, sC_P, img_runStart_C] = CalcStretchStress(rawData, 3);
img_runStart(3) = img_runStart_C;

close(figure(7))

% pad ends of all data arrays to make them the same length (for concatenation):
maxL = max([length(sA_P),length(sB_P),length(sC_P)]);
sA_P(end+1:maxL) = 0;
sA_lam(end+1:maxL) = 0;
sB_P(end+1:maxL) = 0;
sB_lam(end+1:maxL) = 0;
sC_P(end+1:maxL) = 0;
sC_lam(end+1:maxL) = 0;

% output final processed data in cell array (calculated stresses and stretches):
calcData = [{[sA ' stretch'],[sA ' stress (kPa)'],[sB ' stretch'],[sB ' stress (kPa)'],...
    [sC ' stretch'],[sC ' stress (kPa)']};num2cell(sA_lam),num2cell(sA_P),...
    num2cell(sB_lam),num2cell(sB_P),num2cell(sC_lam),num2cell(sC_P)];

% output useful run and measurement information as a cell array:
runInfo = [{'A','B','C','thickness','k_b (mN/mm)',...
    'First image (#)','Run start (A)','Run start (B)','Run start (C)'};...
    {sA,sB,sC,thickness,k_b,firstImg,img_runStart(1),img_runStart(2),img_runStart(3)}]';

% update message:
disp('Analysis complete! Remember to save data.');






%% Internal functions:

function [point_cur] = ShiftROI_ssim(sampleIndex, measureIndex, imageIndex, ROIradius)
    % load previous and current images:
    img_prev = imread(fullfile(pathName,fileNames{firstImg+imageIndex-2}));
    img_cur = imread(fullfile(pathName,fileNames{firstImg+imageIndex-1}));
    % isolate 3D ROI around previous data point in both images:
    if measureIndex == 1
        ROIradius = ROIradius * 3;
        channels = 1;
    else
        channels = [1 2 3];
    end
    ROI_prev = img_prev(M(sampleIndex,measureIndex,imageIndex-1,2)-ROIradius:M(sampleIndex,measureIndex,imageIndex-1,2)+ROIradius, ...
            M(sampleIndex,measureIndex,imageIndex-1,1)-ROIradius:M(sampleIndex,measureIndex,imageIndex-1,1)+ROIradius, channels);
    ROI_cur = img_cur(M(sampleIndex,measureIndex,imageIndex-1,2)-ROIradius:M(sampleIndex,measureIndex,imageIndex-1,2)+ROIradius, ...
        M(sampleIndex,measureIndex,imageIndex-1,1)-ROIradius:M(sampleIndex,measureIndex,imageIndex-1,1)+ROIradius, channels);
    % calculate SSIM between previous and current ROIs:
    ssim_cur = ssim(ROI_cur, ROI_prev);
    % start by assuming this SSIM is the maximum (i.e. closest to 1):
    ssim_max = ssim_cur;
    point_cur = squeeze(M(sampleIndex,measureIndex,imageIndex-1,:));
    % shift ROI vertically to find position of max SSIM:
    for i_pos = 1:length_ROItranslation
        % isolate new shifted 3D ROI in current image:
        ROI_cur = img_cur(M(sampleIndex,measureIndex,imageIndex-1,2)-ROIradius-i_pos:M(sampleIndex,measureIndex,imageIndex-1,2)+ROIradius-i_pos, ...
            M(sampleIndex,measureIndex,imageIndex-1,1)-ROIradius:M(sampleIndex,measureIndex,imageIndex-1,1)+ROIradius, channels);
        % calculate SSIM between previous and new current ROIs: 
        ssim_cur = ssim(ROI_cur, ROI_prev);
        if ssim_cur >= ssim_max
            ssim_max = ssim_cur;
            point_cur = squeeze(M(sampleIndex,measureIndex,imageIndex-1,:)) - [0; i_pos];
        end
    end
end

function [stretch, stress, runStartImage] = CalcStretchStress(data, sample)
    % 1. Determine all the lengths (D & l) from raw data (mm units):
    displacement = squeeze(data(sample,2,img_firstMotion:end)-data(1,1,img_firstMotion:end)) - squeeze(data(sample,2,img_firstMotion)-data(1,1,img_firstMotion)); % calculate displacement between actuator and cantilever
    sampleLength = squeeze(abs(data(sample,3,img_firstMotion:end)-data(sample,4,img_firstMotion:end))); % tissue sample length
    % 2. Set increasing force constraint, then set new initial displacement to 0:
    for i_image = 1:length(displacement)
        fit = fitlm(linspace(0,1,valRange)', displacement(i_image:i_image+valRange-1), 'poly1'); % linear regression
        if fit.Coefficients.Estimate(2) > 0 && fit.Coefficients.pValue(2) < 0.05 % i.e. the slope is confidently (95%) positive in a window of "valRange" number of measurements
            % set measurement i_image as the first measurement (when displacement = 0):
            displacement = squeeze(data(sample,2,img_firstMotion+i_image-1:end)-data(1,1,img_firstMotion+i_image-1:end)) - squeeze(data(sample,2,img_firstMotion+i_image-1)-data(1,1,img_firstMotion+i_image-1));
            sampleLength = squeeze(abs(data(sample,3,img_firstMotion+i_image-1:end)-data(sample,4,img_firstMotion+i_image-1:end)));
            % open image i_image to select initial width for this sample:
            image = imread(fullfile(pathName,fileNames{firstImg+i_image-1}));
            figure(7)
            imshow(image); hold on
            title('Select one point on each side of sample A (left) to measure width.');
            for i_widthPoint = 1:2 % manually select points
                [x,y] = ginput(1);
                width_temp(i_widthPoint) = x;
                scatter(x,y,60,'r','+');
            end
            pause(0.2);
            delete(gca);
            runStartImage = firstImg+i_image-1; % output the starting image
            break
        end
    end
    width = abs(width_temp(1)-width_temp(2)) * res; % initial sample width
    % 3. Calculate force, than calculate stress and stretch:
    force = k_b*displacement; % force on cantilever end
    stress = force / (thickness*width); % P, engineering stress (kPa)
    stretch = sampleLength/sampleLength(1); % λ, stretch (unitless)
end



end

