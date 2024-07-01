function [calcData,rawData,runInfo] = T4Axio_tube_3point()
%% Inputs
% Start by filling in these parameters, then run the function.

% I. experimental setup:
R = 0.064; % cantilever radius (units: mm)
L = 15; % cantilever length (units: mm)
v = 0.019; % actuator velocity (units: mm/s)
res = 0.0114; % image resolution (units: mm/px) (e.g. 0.0091mm/px for 5X zoom)
imagingInterval = 15; % duration between each image in the timelapse (units: s)
rad_i = 0.2; % get the (mean) thickness (t) of samples as a number (mm)
rad_o = 0.4;

% II. analysis settings:
img_start = [8 6 7]; % first img index for each sample at which tube is straight
thldR = 0.01; % threshold level for red channel (0 to 1, higher is more stringent)
thldG = 0.03; % threshold level for green channel (0 to 1, higher is more stringent)
pixelRange_R = [2000,150000]; % [min,max] number of pixels allowed in a red blob (cantilevers, base)
pixelRange_G = [250,5000]; % [min,max] number of pixels allowed in a green blob (tissue markers)


%% Setup

% I. allow the user to select an image file (start of run, usually image 3):
[fileName, pathName] = uigetfile('*.tif', 'Open first image of run:');
fullName = fullfile(pathName,fileName); % these four lines are from Meshoo on MATLAB Answers...
filelist = dir([fileparts(fullName) filesep '*.tif']);
fileNames = {filelist.name}';
imagesAll = numel(filelist);
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

% IV. set stringency of positive slope test for finding start of run:
valRange = 5; % X number of measurements in a row to be tested for (+) slope

% V. initialize a 3D array (m) to hold measured values (and an array (zones) to hold zone bounds):
m = zeros(3,5,images/3); % dimensions (in order): sample (A,B,C), measurement, image number

% VI. use Inputs to calculate bending stiffness and cross-sectional area:
E = 410e9; % known modulus (E) of tungsten (Pa)
k_b = (3*E*pi*(R/1000)^4)/(4*(L/1000)^3);
csArea = 2*pi*(rad_o^2-rad_i^2); % 2x because there are two side of tube

% VII. welcome message:
disp('Welcome!');

% some variables that affect the display of example images while the analysis is happening:
checkImageEvery = 4;


%% Image Processing:

figure(1) % to display blob-analysis example images as processing progresses
figure(2) % thresholded red images
figure(3) % thresholded green images
figure(4) % to display point measurement example images as processing progresses
interationCount = 0;
for i = 1:3:images
    interationCount = interationCount + 1;
    
    % I. import the next image into MATLAB:
    imgNum = firstImg + i - 1;
    imgR = imread(fullfile(pathName,fileNames{imgNum}));
    imgG = imread(fullfile(pathName,fileNames{imgNum+1}));
    img = imread(fullfile(pathName,fileNames{imgNum+2}));
    
    % II. split RGB channels, theshold:
    imgRthr = imbinarize(imgR,thldR);
    imgGthr = imbinarize(imgG,thldG);
    
    % III. identify and measure blobs:
    % eliminate "small" and "large" blobs:
    imgRfiltered = bwareafilt(imgRthr,pixelRange_R);
    imgGfiltered = bwareafilt(imgGthr,pixelRange_G);
    % use the largest blob (cantilever "block") to black out upper part of image:
    if i == 1
        figure(5)
        imshow(img+50); hold on
        title('Select a point just above the cantilvers but below the actuator head:');
        [~,pointBlock] = ginput(1);
        close(5);
    else
        pointBlock = pointBlock - v * imagingInterval / res;
        if pointBlock < 1
            pointBlock = 1;
        end
    end
    imgRfiltered(1:round(pointBlock),:) = 0; % black out all pixels in that area of image, using centroid from R or G image (whichever is more stringent)
    imgGfiltered(1:round(pointBlock),:) = 0;
    % label and count blobs:
    [labelsR,numR] = bwlabel(imgRfiltered);
    [labelsG,numG] = bwlabel(imgGfiltered);
    % calculate centroids to get data points (then convert struct to matrix form):
    pointsR_struct = regionprops(labelsR,'Centroid');
    pointsG_struct = regionprops(labelsG,'Centroid');
    pointsR = zeros(numR,2); % initialize variables...
    pointsG = zeros(numG,2);
    for j = 1:numR
        pointsR(j,1:2) = [pointsR_struct(j).Centroid(1),pointsR_struct(j).Centroid(2)];
    end
    for j = 1:numG
        pointsG(j,1:2) = [pointsG_struct(j).Centroid(1),pointsG_struct(j).Centroid(2)];
    end
    
    % IV. displaying blob-analysis example images, show every 4rth image:
    if rem(imgNum-firstImg,checkImageEvery) == 0 && (imgNum-firstImg)/checkImageEvery+1 <= 5*6
        figure(1)
        subplot(5,6,(imgNum-firstImg)/checkImageEvery+1);
        imshow(img);
        figure(2)
        subplot(5,6,(imgNum-firstImg)/checkImageEvery+1);
        imshow(imgRfiltered); hold on
        figure(3)
        subplot(5,6,(imgNum-firstImg)/checkImageEvery+1);
        imshow(imgGfiltered); hold on
    end
    
    % V. average centroids of nearest neighbor blobs until number "red" point = 4 and number "green" points = 9:
    % red points in the bottom 1/3 should be cantilever:
    [indices] = find(pointsR(:,2) > length(imgRfiltered) * (2/3)); % points in the bottom 1/3 of image
    xNew = mean(pointsR(indices(:),1));
    yNew = mean(pointsR(indices(:),2));
    pointsR(indices(:),:) = []; % remove these points from the points list
    pointsR(end+1,:) = [xNew,yNew]; % add the new averaged point to the points list as a centroid with new x and y
    pointsR = pointsR(~isnan(pointsR(:,1)),:); % remove nans, just in case
    % if run had actuator tracker attached, uncomment this to remove left-most red blob (for consistency):
    [~,index] = min(pointsR(:,1));
    pointsR(index,:) = [];
    % rest of red points:
    while length(pointsR) > 4
        nearestPnts = [0,0]; % [point 1 index, point 2 index]
        leastDist = length(img)^2; % initialize variable to hold smallest measured distance (defined here in such a way as to ensure that it always starts bigger than any possible distance within the image)
        for j = 1:length(pointsR)
            for k = 1:length(pointsR)
                if j ~= k
                    dist = (pointsR(j,1) - pointsR(k,1))^2 +...
                        (pointsR(j,2) - pointsR(k,2))^2; % sqrt() skipped to save time
                    if dist < leastDist
                        leastDist = dist;
                        nearestPnts = [j,k];
                    end
                end
            end
        end
        xNew = mean([pointsR(nearestPnts(1),1),pointsR(nearestPnts(2),1)]); % calculate average x position
        yNew = mean([pointsR(nearestPnts(1),2),pointsR(nearestPnts(2),2)]); % calculate average y position
        pointsR([nearestPnts(1),nearestPnts(2)],:) = []; % remove the two nearest points from the points list
        pointsR(end+1,:) = [xNew,yNew]; % add the new averaged point to the points list as a centroid with new x and y
    end
    while length(pointsG) > 9
        nearestPnts = [0,0]; % point 1 index, point 2 index
        leastDist = length(img)^2; % initialize variable to hold smallest measured distance
        for j = 1:length(pointsG)
            for k = 1:length(pointsG)
                if j ~= k
                    dist = (pointsG(j,1) - pointsG(k,1))^2 +...
                        (pointsG(j,2) - pointsG(k,2))^2; % sqrt() skipped to save time
                    if dist < leastDist
                        leastDist = dist;
                        nearestPnts = [j,k];
                    end
                end
            end
        end
        xNew = mean([pointsG(nearestPnts(1),1),pointsG(nearestPnts(2),1)]);
        yNew = mean([pointsG(nearestPnts(1),2),pointsG(nearestPnts(2),2)]);
        pointsG([nearestPnts(1),nearestPnts(2)],:) = [];
        pointsG(end+1,:) = [xNew,yNew];
    end
    if length(pointsR) < 4 || length(pointsG) < 9
        disp(['Processing ended early due to too few points in image. Got to image ' num2str(interationCount) ' of ' num2str(images/3) '. Points (R,G): ' num2str(length(pointsR)) ', ' num2str(length(pointsG))]);
        break
    end
    
    % VI. for visual confirmation, show every 4th image with measured points:
    % display raw images with collected data points (in figure(4)), just to
    % confirm things are working correctly (see subsection IX in mT4auto.m
    % Image Processing)
%     if rem(imgNum-firstImg,checkImageEvery) == 0
%         figure(4)
%         subplot(5,6,(imgNum-firstImg)/checkImageEvery+1);
%         imshow(img); hold on
%         for j = 1:length(pointsR)
%             scatter(pointsR(j,1),pointsR(j,2),250,'c','+');
%             scatter(pointsR(j,1),pointsR(j,2),150,'c','.');
%         end
%         for j = 1:length(pointsG)
%             scatter(pointsG(j,1),pointsG(j,2),250,'y','+');
%             scatter(pointsG(j,1),pointsG(j,2),150,'y','.');
%         end
%     end
    
    if rem(imgNum-firstImg,checkImageEvery) == 0 && (imgNum-firstImg)/checkImageEvery+1 <= 5*6
        figure(4)
        subplot(5,6,(imgNum-firstImg)/checkImageEvery+1);
        imshow(img); hold on
        
        % VII. for Red points, use relative point positions to associate points with corresponding structures:
        points = zeros(3,5); % for each sample (1-3,:), holds y positions of red points (base(i,1) and cantilever (i,2)) and green points (i,3:5)
        % cantilever (lowest point):
        [~,index] = max(pointsR(:,2));
        scatter(pointsR(index,1),pointsR(index,2),100,'w','+');
        points(1,2) = pointsR(index,2);
        points(2,2) = pointsR(index,2);
        points(3,2) = pointsR(index,2);
        pointsR(index,:) = []; % remove that point from the pointsR array
        % other points:
        pointsR = sortrows(pointsR);
        points(1,1) = pointsR(1,2);
        points(2,1) = pointsR(2,2);
        points(3,1) = pointsR(3,2);
        scatter(pointsR(1,1),pointsR(1,2),50,'w','x');
        scatter(pointsR(2,1),pointsR(2,2),50,'y','x');
        scatter(pointsR(3,1),pointsR(3,2),50,'c','x');

        % VIII. for Green points, use relative point positions to associate points with corresponding structures (^ see VII for line-by-line details):
        % for the 12 green points, use x values to identify which points go to which sample and which side of which sample (save the paired y values as data):
        % (Note: top vs bottom point (from y value) it not important here b/c we will only use these points to measure the absolute y distance between them for strain calculations;
        % however, we do want to keep track of the left vs right set of points for each sample so that we can calculate the strain on the left and right part of tube separately)
        pointsG = sortrows(pointsG);
        points(1,3:5) = pointsG(1:3,2);
        points(2,3:5) = pointsG(4:6,2);
        points(3,3:5) = pointsG(7:9,2);
        scatter(pointsG(1,1),pointsG(1,2),50,'w','<');
        scatter(pointsG(2,1),pointsG(2,2),50,'w','<');
        scatter(pointsG(3,1),pointsG(3,2),50,'w','<');
        scatter(pointsG(4,1),pointsG(4,2),50,'y','<');
        scatter(pointsG(5,1),pointsG(5,2),50,'y','<');
        scatter(pointsG(6,1),pointsG(6,2),50,'y','<');
        scatter(pointsG(7,1),pointsG(7,2),50,'c','<');
        scatter(pointsG(8,1),pointsG(8,2),50,'c','<');
        scatter(pointsG(9,1),pointsG(9,2),50,'c','<');
    else
        % VII. for Red points, use relative point positions to associate points with corresponding structures:
        points = zeros(3,4); % for each sample (1-3,:), holds y positions of red points (base(i,1) and cantilever (i,2)) and green points (left (i,3-4) and right (i,5-6))
        % cantilever (lowest point):
        [~,index] = max(pointsR(:,2));
        points(1,2) = pointsR(index,2);
        points(2,2) = pointsR(index,2);
        points(3,2) = pointsR(index,2);
        pointsR(index,:) = []; % remove that point from the pointsR array
        % other points:
        pointsR = sortrows(pointsR);
        points(1,1) = pointsR(1,2);
        points(2,1) = pointsR(2,2);
        points(3,1) = pointsR(3,2);

        % VIII. for Green points, use relative point positions to associate points with corresponding structures (^ see VII for line-by-line details):
        % for the 12 green points, use x values to identify which points go to which sample and which side of which sample (save the paired y values as data):
        % (Note: top vs bottom point (from y value) it not important here b/c we will only use these points to measure the absolute y distance between them for strain calculations;
        % however, we do want to keep track of the left vs right set of points for each sample so that we can calculate the strain on the left and right part of tube separately)
        % left green pair for sample A:
        pointsG = sortrows(pointsG);
        points(1,3:5) = pointsG(1:3,2);
        points(2,3:5) = pointsG(4:6,2);
        points(3,3:5) = pointsG(7:9,2);
    end
    
    % IX. save values to m as distance measurments (mm units):
    m(:,:,interationCount) = points * res;
end

% X. output raw data (distance measurements only):
m(:,3:5,:) = sort(m(:,3:5,:), 2);
rawData = m;

% XI. update message:
disp('Measurements complete. Calculations underway...');

%% Calculations:

% I. define EXPECTED cantilever displacement (based on v):
time = zeros(images/3,1); % define number of timepoints based on number of images
for i = 1:images/3-1
    time(i+1) = time(i)+imagingInterval; % set current time (from t=0) at each timepoint (15sec between images)
end
D = v*time; % (displacement) = (run velocity) * (current time)

% II. Cantilever A:
% i. Determine all the lengths (d, D, l) from raw data (mm units):
sA_d = squeeze(abs(m(1,1,:)-m(1,2,:))); % cantilever distance from base (squeeze() removes singleton dimensions)
sA_D = (sA_d(1)+D)-sA_d; % calculated displacement
sA_l_upper = squeeze(abs(m(1,3,:)-m(1,4,:))); % sample length, left tube
sA_l_lower = squeeze(abs(m(1,4,:)-m(1,5,:))); % sample length, left tube
sA_l = (sA_l_upper + sA_l_lower) / 2; % mean length
% ii. apply 'start image' constraint:
sA_d = sA_d(img_start(1):end);
sA_D = (sA_d(1)+D(1:end+1-img_start(1)))-sA_d;
sA_l = sA_l(img_start(1):end);
% iii. Set increasing force constraint, then set new initial displacement to 0:
for i = 1:length(sA_D)
    p = polyfit(linspace(0,1,valRange)',sA_D(i:i+valRange-1),1); % linear regression, p(1) = slope
    if p(1) > 0 % i.e. the regression slope is positive in a window of "valRange" number of measurements
        % set measurement i as the first measurement (when displacement = 0):
        sA_d = sA_d(i:end);
        sA_D = (sA_d(1)+D(1:end+2-img_start(1)-i))-sA_d;
        sA_l = sA_l(i:end);
        break
    end
end
% iv. Calculate force, than calculate stress and strain:
sA_F = k_b*sA_D; % force on cantilever end
sA_s = sA_F/csArea; % s, engineering stress (kPa)
sA_e = (sA_l-sA_l(1))/sA_l(1); % e, engineering strain (unitless)

% III. Cantilever B:
% I. Determine all the lengths (d, D, l) from raw data (mm units):
sB_d = squeeze(abs(m(2,1,:)-m(2,2,:))); % cantilever distance from base (squeeze() removes singleton dimensions)
sB_D = (sB_d(1)+D)-sB_d; % calculated displacement
sB_l_upper = squeeze(abs(m(2,3,:)-m(2,4,:))); % sample length, left tube
sB_l_lower = squeeze(abs(m(2,4,:)-m(2,5,:))); % sample length, left tube
sB_l = (sB_l_upper + sB_l_lower) / 2; % mean length
% ii. apply 'start image' constraint:
sB_d = sB_d(img_start(2):end);
sB_D = (sB_d(1)+D(1:end+1-img_start(2)))-sB_d;
sB_l = sB_l(img_start(2):end);
% iii. Set increasing force constraint, then set new initial displacement to 0:
for i = 1:length(sB_D)
    p = polyfit(linspace(0,1,valRange)',sB_D(i:i+valRange-1),1); % linear regression, p(1) = slope
    if p(1) > 0 % i.e. the regression slope is positive in a window of "valRange" number of measurements
        % set measurement i as the first measurement (when displacement = 0):
        sB_d = sB_d(i:end);
        sB_D = (sB_d(1)+D(1:end+2-img_start(2)-i))-sB_d;
        sB_l = sB_l(i:end);
        break
    end
end
% iv. Calculate force, than calculate stress and strain:
sB_F = k_b*sB_D; % force on cantilever end
sB_s = sB_F/csArea; % s, engineering stress (kPa)
sB_e = (sB_l-sB_l(1))/sB_l(1); % e, engineering strain (unitless)

% IV. Cantilever C:
sC_d = squeeze(abs(m(3,1,:)-m(3,2,:))); % cantilever distance from base (squeeze() removes singleton dimensions)
sC_D = (sC_d(1)+D)-sC_d; % calculated displacement
sC_l_upper = squeeze(abs(m(3,3,:)-m(3,4,:))); % sample length, left tube
sC_l_lower = squeeze(abs(m(3,4,:)-m(3,5,:))); % sample length, left tube
sC_l = (sC_l_upper + sC_l_lower) / 2; % mean length
% ii. apply 'start image' constraint:
sC_d = sC_d(img_start(3):end);
sC_D = (sC_d(1)+D(1:end+1-img_start(3)))-sC_d;
sC_l = sC_l(img_start(3):end);
% iii. Set increasing force constraint, then set new initial displacement to 0:
for i = 1:length(sC_D)
    p = polyfit(linspace(0,1,valRange)',sC_D(i:i+valRange-1),1); % linear regression, p(1) = slope
    if p(1) > 0 % i.e. the regression slope is positive in a window of "valRange" number of measurements
        % set measurement i as the first measurement (when displacement = 0):
        sC_d = sC_d(i:end);
        sC_D = (sC_d(1)+D(1:end+2-img_start(3)-i))-sC_d;
        sC_l = sC_l(i:end);
        break
    end
end
% iv. Calculate force, than calculate stress and strain:
sC_F = k_b*sC_D; % force on cantilever end
sC_s = sC_F/csArea; % s, engineering stress (kPa)
sC_e = (sC_l-sC_l(1))/sC_l(1); % e, engineering strain (unitless)

% V. pad ends of all data arrays to make them the same length (for concatenation):
maxL = max([length(sA_s),length(sB_s),length(sC_s)]);
sA_s(end+1:maxL) = 0;
sA_e(end+1:maxL) = 0;
sB_s(end+1:maxL) = 0;
sB_e(end+1:maxL) = 0;
sC_s(end+1:maxL) = 0;
sC_e(end+1:maxL) = 0;

% VI. output final processed data in cell array (calculated stresses and strains):
calcData = [{[sA ' strain'],[sA ' stress (kPa)'],[sB ' strain'],[sB ' stress (kPa)'],...
    [sC ' strain'],[sC ' stress (kPa)']};num2cell(sA_e),num2cell(sA_s),...
    num2cell(sB_e),num2cell(sB_s),num2cell(sC_e),num2cell(sC_s)];

% VII. output useful run and measurement information as a cell array:
runInfo = [{'A','B','C','inner r (mm)','outer r (mm)','k_b (mN/mm)','v (mm/s)',...
    'First image (#)','Meas. interval (s)','Threshold (R)','Threshold (G)','Bolb size range (R)','Bolb size range (G)'};...
    {sA,sB,sC,rad_i,rad_o,k_b,v,firstImg,imagingInterval,thldR,thldG,pixelRange_R,pixelRange_G}]';

% VIII. update message:
disp('Analysis complete! Remember to save data.');


