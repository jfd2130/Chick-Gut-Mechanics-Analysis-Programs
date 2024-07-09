function [calcData,rawData,runInfo] = MT3Axio_mesentaryAuto()
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

% IV. set image analysis parameters [MANUALLY]:
thldR = 0.1; % threshold level (0 to 1, higher is more stringent)
thldG = 0.035;
pixelRange_R = [500,150000]; % min and max number of pixels in a red blob (cantilevers, base)
pixelRange_G = [1000,10000]; % min and max number of pixels in a green blob (tissue markers)

% V. set stringency of positive slope test for finding start of run:
valRange = 5; % X number of measurements in a row to be tested for (+) slope

% VI. initialize a 3D array (m) to hold measured values (and an array (zones) to hold zone bounds):
m = zeros(3,4,images/3); % dimensions (in order): sample (A,B,C), measurement, image number

% VII. welcome message:
disp('Welcome!');

checkImageEvery = 4;
brightenImgBy = 50;

%% Initialization:

% I. present dialogue box to enter cantilever and velocity parameters for the run:
parameters = inputdlg({'Cantilever radius (mm):','Cantilver length (mm):',...
    'Velocity (mm/s):'},'Set mechanical parameters:',[1 60]);
R = str2double(parameters{1})*10^-3; % get the radius (R) as a number (m)
L = str2double(parameters{2})*10^-3; % get the length (L) as a number (m)
v = str2double(parameters{3}); % get the velocity (v) as a number (mm/s)
E = 410e9; % known modulus (E) of tungsten (Pa)
% use R, L, and E to calculate bending stiffness:
k_b = (3*E*pi*R^4)/(4*L^3);

% II. get the resolution:
resolution = inputdlg({'resolution (mm/px):'},'Set image resolution:',[1 60]);
res = str2double(resolution{1});
imagingInterval = 20; % seconds between images

% III. present dialogue box to enter sample radii:
thicknessdlg = inputdlg({'thickness'},'Set samples thickness:',[1 60]);
thickness = str2double(thicknessdlg{1}); % get the thickness (t) of samples on A as a number (mm)

% IV. update message:
disp('Initialization complete. Image processing underway...');


%% Image Processing

figure(1) % to display blob-analysis example images as processing progresses
figure(2) % thresholded red images
figure(3) % thresholded green images
figure(4) % to display point measurement example images as processing progresses
for i = 1:images/3
    imgNum = i*3-2;
    % I. import images in 3 channels 
    imgBf = imread(fullfile(pathName,fileNames{imgNum}));
    imgG = imread(fullfile(pathName,fileNames{imgNum+1}));
    imgR = imread(fullfile(pathName,fileNames{imgNum+2}));
    
    % II. split RGB channels, theshold:
    imgRthr = imbinarize(imgR,thldR);
    imgGthr = imbinarize(imgG,thldG);
    
    % III. identify and measure blobs:
    % eliminate "small" and "large" blobs:
    imgRfiltered = bwareafilt(imgRthr,pixelRange_R);
    imgGfiltered = bwareafilt(imgGthr,pixelRange_G);
    % use the largest blob (cantilever "block") to black out upper part of image:
    if i == 1
        imgRblock = bwareafilt(imgRthr,1); % keep only the one largest blob
        [labelRblock,~] = bwlabel(imgRblock);
        pointRblock = regionprops(labelRblock,'Centroid');
        imgGblock = bwareafilt(imgGthr,1); % keep only the one largest blob
        [labelGblock,~] = bwlabel(imgGblock);
        pointGblock = regionprops(labelGblock,'Centroid');
        if pointRblock(1).Centroid(2) >= pointGblock(1).Centroid(2)
            pointBlock = pointRblock(1).Centroid(2) * 2;
        else
            pointBlock = pointGblock(1).Centroid(2) * 2;
        end
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
    
    % IV. displaying blob-analysis example images, show every 4th image:
    if rem(i,checkImageEvery) == 0 && (i)/checkImageEvery+1 <= 5*6
        figure(1)
        subplot(5,6,(i)/checkImageEvery);
        imshow(imgBf+brightenImgBy);
        figure(2)
        subplot(5,6,(i)/checkImageEvery);
        imshow(imgRfiltered); hold on
        scatter(round(length(imgRfiltered)/2),pointBlock,250,'r','x')
        figure(3)
        subplot(5,6,(i)/checkImageEvery);
        imshow(imgGfiltered); hold on
        scatter(round(length(imgGfiltered)/2),pointBlock,250,'g','x')
    end
    
    
    % V. average centroids of nearest neighbor blobs until number "red" point = 4 and number "green" points = 6:
    while length(pointsR) > 4
        nearestPnts = [0,0]; % [point 1 index, point 2 index]
        leastDist = length(imgBf)^2; % initialize variable to hold smallest measured distance (defined here in such a way as to ensure that it always starts bigger than any possible distance within the image)
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
    while length(pointsG) > 6
        nearestPnts = [0,0]; % point 1 index, point 2 index
        leastDist = length(imgBf)^2; % initialize variable to hold smallest measured distance
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
    if length(pointsR) < 4 || length(pointsG) < 6
        disp(['Processing ended early due to too few points in image. Got to image ' num2str(i) ' of ' num2str(images/3) '. Points (R,G): ' num2str(length(pointsR)) ', ' num2str(length(pointsG))]);
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
    
    if rem(i,checkImageEvery) == 0 && (i)/checkImageEvery+1 <= 5*6
        figure(4)
        subplot(5,6,(i)/checkImageEvery);
        imshow(imgBf+brightenImgBy); hold on
        
        % VII. for Red points, use relative point positions to associate points with corresponding structures:
        points = zeros(3,4); % for each sample (1-3,:), holds y positions of red points (base(i,1) and cantilever (i,2)) and green points (left (i,3) and right (i,4))
        %find cantilever
        [~,index] = max(pointsR(:,2));
        points(:,1) = pointsR(index,2);
        scatter(pointsR(index,1),pointsR(index,2),50,'w','+');
        pointsR(index,:) = [];
        % red pair for sample A:
        [~,index] = min(pointsR(:,1)); % get index of current farthest left red point, this belongs to the left sample (1)
        points(1,2) = pointsR(index,2); % save the y value of that point
        scatter(pointsR(index,1),pointsR(index,2),50,'w','x');
        pointsR(index,:) = []; % remove that point from the pointsR array 
        % red pair for sample B:
        [~,index] = min(pointsR(:,1)); % get index of current farthest left red point, this belongs to the center sample (2)
        points(2,2) = pointsR(index,2); % save the y value of that point
        scatter(pointsR(index,1),pointsR(index,2),50,'y','x');
        pointsR(index,:) = []; % remove that point from the pointsR array
        % red pair for sample C:
        scatter(pointsR(1,1),pointsR(1,2),50,'c','x');
        points(3,2) = pointsR(end,2); % save the y values of the remaining two points, these belong to the right sample (C)

        % VIII. for Green points, use relative point positions to associate points with corresponding structures (^ see VII for line-by-line details):
        % for the 6 green points, use x values to identify which points go to which sample and which side of which sample (save the paired y values as data):
        % (Note: top vs bottom point (from y value) it not important here b/c we will only use these points to measure the absolute y distance between them for strain calculations;
        % however, we do want to keep track of the left vs right set of points for each sample so that we can calculate the strain on the left and right part of tube separately)
        % left green point for sample A:
        [~,index] = min(pointsG(:,1));
        points(1,3) = pointsG(index,2);
        scatter(pointsG(index,1),pointsG(index,2),50,'w','<');
        pointsG(index,:) = [];
        [~,index] = min(pointsG(:,1));
        points(1,4) = pointsG(index,2);
        scatter(pointsG(index,1),pointsG(index,2),50,'w','<');
        pointsG(index,:) = [];

        % left green point for sample B:
        [~,index] = min(pointsG(:,1));
        points(2,3) = pointsG(index,2);
        scatter(pointsG(index,1),pointsG(index,2),50,'y','<');
        pointsG(index,:) = [];
        [~,index] = min(pointsG(:,1));
        points(2,4) = pointsG(index,2);
        scatter(pointsG(index,1),pointsG(index,2),50,'y','<');
        pointsG(index,:) = [];

        % left green point for sample C:
        scatter(pointsG(1,1),pointsG(1,2),50,'c','<');
        scatter(pointsG(2,1),pointsG(2,2),50,'c','<');
        points(3,3:4) = pointsG(:,2)';
    else
        % VII. for Red points, use relative point positions to associate points with corresponding structures:
        points = zeros(3,4); % for each sample (1-3,:), holds y positions of red points (base(i,1) and cantilever (i,2)) and green points (left (i,3) and right (i,4))
        %find cantilever
        [~,index] = max(pointsR(:,2));
        points(:,1) = pointsR(index,2);
        pointsR(index,:) = [];
        % red pair for sample A:
        [~,index] = min(pointsR(:,1)); % get index of current farthest left red point, this belongs to the left sample (1)
        points(1,2) = pointsR(index,2); % save the y value of that point
        pointsR(index,:) = []; % remove that point from the pointsR array 
        % red pair for sample B:
        [~,index] = min(pointsR(:,1)); % get index of current farthest left red point, this belongs to the center sample (2)
        points(2,2) = pointsR(index,2); % save the y value of that point
        pointsR(index,:) = []; % remove that point from the pointsR array
        % red pair for sample C:
        points(3,2) = pointsR(end,2); % save the y values of the remaining two points, these belong to the right sample (C)

        % VIII. for Green points, use relative point positions to associate points with corresponding structures (^ see VII for line-by-line details):
        % for the 6 green points, use x values to identify which points go to which sample and which side of which sample (save the paired y values as data):
        % (Note: top vs bottom point (from y value) it not important here b/c we will only use these points to measure the absolute y distance between them for strain calculations;
        % however, we do want to keep track of the left vs right set of points for each sample so that we can calculate the strain on the left and right part of tube separately)
        % left green point for sample A:
        [~,index] = min(pointsG(:,1));
        points(1,3) = pointsG(index,2);
        pointsG(index,:) = [];
        [~,index] = min(pointsG(:,1));
        points(1,4) = pointsG(index,2);
        pointsG(index,:) = [];

        % left green point for sample B:
        [~,index] = min(pointsG(:,1));
        points(2,3) = pointsG(index,2);
        pointsG(index,:) = [];
        [~,index] = min(pointsG(:,1));
        points(2,4) = pointsG(index,2);
        pointsG(index,:) = [];

        % left green point for sample C:
        points(3,3:4) = pointsG(:,2)';
    end
    
    % IX. save values to m as distance measurments (mm units):
    m(:,:,i) = points * res;
end

% X. output raw data (distance measurements only):
rawData = m;

% XI. update message:
disp('Measurements complete. Calculations underway...');

%% Calculations:

figure(5)

% I. define EXPECTED cantilever displacement (based on v):
time = zeros(images/3,1); % define number of timepoints based on number of images
for i = 1:images/3-1
    time(i+1) = time(i)+20; % set current time (from t=0) at each timepoint (20sec between images)
end
D = v*time; % (displacement) = (run velocity) * (current time)

% II. Cantilever A:
% i. Determine all the lengths (d, D, l) from raw data (mm units):
sA_d = squeeze(abs(m(1,1,:)-m(1,2,:))); % cantilever distance from base (squeeze() removes singleton dimensions)
sA_D = (sA_d(1)+D)-sA_d; % calculated displacement
sA_l = squeeze(abs(m(1,3,:)-m(1,4,:))); % sample length
% ii. Set increasing force constraint, then set new initial displacement to 0:
for i = 1:length(sA_D)
    p = polyfit(linspace(0,1,valRange)',sA_D(i:i+valRange-1),1); % linear regression, p(1) = slope
    if p(1) > 0 % i.e. the regression slope is positive in a window of "valRange" number of measurements
        % set measurement i as the first measurement (when displacement = 0):
        sA_d = sA_d(i:end);
        sA_D = (sA_d(1)+D(1:end+1-i))-sA_d;
        sA_l = sA_l(i:end);
        % open image i to select initial width for this sample:
        img = imread(fullfile(pathName,fileNames{firstImg+i*3-3}));
        figure(5)
        imshow(img+brightenImgBy); hold on
        title('Select one point on each side of sample A (left) to measure width.');
        for j = 1:2 % manually select points
            [x,y] = ginput(1);
            wTemp(j) = x;
            scatter(x,y,60,'r','+');
        end
        delete(gca);
        break
    end
end
sA_w = abs(wTemp(1)-wTemp(2))*res; % initial sample width
% iii. Calculate force, than calculate stress and strain:
sA_F = k_b*sA_D; % force on cantilever end
sA_s = sA_F/(thickness*sA_w); % s, engineering stress (kPa)
sA_e = (sA_l-sA_l(1))/sA_l(1); % e, engineering strain (unitless)


% III. Cantilever B:
% i. Determine all the lengths (d, D, l) from raw data (mm units):
sB_d = squeeze(abs(m(2,1,:)-m(2,2,:))); % cantilever distance from base (squeeze() removes singleton dimensions)
sB_D = (sB_d(1)+D)-sB_d; % calculated displacement
sB_l = squeeze(abs(m(2,3,:)-m(2,4,:))); % sample length
% ii. Set increasing force constraint, then set new initial displacement to 0:
for i = 1:length(sB_D)
    p = polyfit(linspace(0,1,valRange)',sB_D(i:i+valRange-1),1); % linear regression, p(1) = slope
    if p(1) > 0 % i.e. the regression slope is positive in a window of "valRange" number of measurements
        % set measurement i as the first measurement (when displacement = 0):
        sB_d = sB_d(i:end);
        sB_D = (sB_d(1)+D(1:end+1-i))-sB_d;
        sB_l = sB_l(i:end);
        % open image i to select initial width for this sample:
        img = imread(fullfile(pathName,fileNames{firstImg+i*3-3}));
        figure(5)
        imshow(img+brightenImgBy); hold on
        title('Select one point on each side of sample A (left) to measure width.');
        for j = 1:2 % manually select points
            [x,y] = ginput(1);
            wTemp(j) = x;
            scatter(x,y,60,'r','+');
        end
        delete(gca);
        break
    end
end
sB_w = abs(wTemp(1)-wTemp(2))*res; % initial sample width
% iii. Calculate force, than calculate stress and strain:
sB_F = k_b*sB_D; % force on cantilever end
sB_s = sB_F/(thickness*sB_w); % s, engineering stress (kPa)
sB_e = (sB_l-sB_l(1))/sB_l(1); % e, engineering strain (unitless)

% IV. Cantilever C:
% i. Determine all the lengths (d, D, l) from raw data (mm units):
sC_d = squeeze(abs(m(3,1,:)-m(3,2,:))); % cantilever distance from base (squeeze() removes singleton dimensions)
sC_D = (sC_d(1)+D)-sC_d; % calculated displacement
sC_l = squeeze(abs(m(3,3,:)-m(3,4,:))); % sample length
% ii. Set increasing force constraint, then set new initial displacement to 0:
for i = 1:length(sC_D)
    p = polyfit(linspace(0,1,valRange)',sC_D(i:i+valRange-1),1); % linear regression, p(1) = slope
    if p(1) > 0 % i.e. the regression slope is positive in a window of "valRange" number of measurements
        % set measurement i as the first measurement (when displacement = 0):
        sC_d = sC_d(i:end);
        sC_D = (sC_d(1)+D(1:end+1-i))-sC_d;
        sC_l = sC_l(i:end);
        % open image i to select initial width for this sample:
        img = imread(fullfile(pathName,fileNames{firstImg+i*3-3}));
        figure(5)
        imshow(img+brightenImgBy); hold on
        title('Select one point on each side of sample A (left) to measure width.');
        for j = 1:2 % manually select points
            [x,y] = ginput(1);
            wTemp(j) = x;
            scatter(x,y,60,'r','+');
        end
        delete(gca);
        close(figure(5));
        break
    end
end
sC_w = abs(wTemp(1)-wTemp(2))*res; % initial sample width
% iii. Calculate force, than calculate stress and strain:
sC_F = k_b*sC_D; % force on cantilever end
sC_s = sC_F/(thickness*sC_w); % s, engineering stress (kPa)
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
runInfo = [{'A','B','C','thickness','k_b (mN/mm)','v (mm/s)',...
    'First image (#)','Meas. interval (s)','Threshold (R)','Threshold (G)','Bolb size range (R)','Bolb size range (G)'};...
    {sA,sB,sC,thickness,k_b,v,firstImg,imagingInterval,thldR,thldG,pixelRange_R,pixelRange_G}]';

% VIII. update message:
disp('Analysis complete! Remember to save data.');


