function [binSums] = fiberAlignment(img_spat)

% References:
% Ayes et al. Measuring fiber alignment in electrospun scaffolds: a user’s guide to the 2D fast Fourier transform approach. J Biomat Sci Polymer Edn, 19(5). 2008.
% Suever & Carlos Borau discussion on Stack Overflow

% allow the user to select an image file for import:
% [fileName, pathName] = uigetfile('*.tif', 'Open image:');
% warning('off','images:initSize:adjustingMag'); % reduces Command Window clutter
% figure

% import the spatial domain image:
%img_spat = imread(fullfile(pathName,fileName)); % import the file
%img_spat = rgb2gray(img_spat); % need gray image
% ensure the image is square:
if length(img_spat(1,:)) > length(img_spat(:,1))
    % if wider than tall, make less wide:
    img_spat = img_spat(:,1:length(img_spat(:,1)));
elseif length(img_spat(1,:)) < length(img_spat(:,1))
    % if taller than wide, make less tall:
    img_spat = img_spat(1:length(img_spat(1,:)),:); 
end
% ensure the image has odd dimesions (so that there is a true center pixel):
if rem(length(img_spat), 2) == 0
    % if dimesions are even, remove a pixel from the width and height:
    img_spat = img_spat(1:end-1,1:end-1);
end
img_spat = double(img_spat);
img_spat = img_spat - min(min(img_spat));
img_spat = img_spat / max(max(img_spat));
% display the image:
figure
subplot(2,3,1); imshow(img_spat,[]); set(gca,'XTick',[],'YTick',[]); title('spatial'); % show image in subplot

% construct polar coordinate map:
polarMap = zeros(length(img_spat),length(img_spat),2); % initialize (stores radius (:,:,1) and angle (:,:,2) from center for each image pixel position)
index_center = (length(img_spat) + 1) / 2; % get position of center pixel (image is square, so x and y coordinates are the same)
for y = 1:length(img_spat)
    for x = 1:length(img_spat)
        % compute and store radius to pixel from center:
        polarMap(y,x,1) = sqrt((y - index_center)^2 + (x - index_center)^2);
        % compute and store angle to pixel from center:
        polarMap(y,x,2) = rad2deg(atan2((y - index_center), (x - index_center)));
    end
end

% apply edge mask to mitigate edge effects:
maskMin = index_center * 0.79; % distance from center that edge gradient to black starts
maskMax = index_center * 0.99; % distance from center beyond which all pixels are made black
img_spat_masked = double(img_spat); % work with doubles (not uint8) when performing mathematical operation on pixel values
for y = 1:length(img_spat_masked)
    for x = 1:length(img_spat_masked)
        if polarMap(y,x,1) > maskMin
            if polarMap(y,x,1) > maskMax
                % make far edges black:
                img_spat_masked(y,x) = 0;
            else
                % create a gradient from black to the image between maskMax and maskMin distances:
                img_spat_masked(y,x) = img_spat_masked(y,x) * (maskMax - polarMap(y,x,1)) / (maskMax - maskMin);
            end
        end
    end
end
% display the masked spatial image:
subplot(2,3,2); imshow(img_spat_masked,[]); set(gca,'XTick',[],'YTick',[]); title('spatial, masked'); % show image in subplot

% convert to frequency domain:
data_freq = fft2(img_spat_masked); % take FFT
data_freq_real = fftshift(abs(data_freq)); % we need to work with the real, shifted version of the frequency domain
data_freq_real = rot90(data_freq_real); % rotate the transform by 90 degrees to realign it with the spatial image
% display the frequency image:
subplot(2,3,3); imshow(data_freq_real,[]); set(gca,'XTick',[],'YTick',[]); title('nat log frequency'); % show image in subplot

% sum frequency values in pixels in degree angle bins:
binSize = 1; % bin size (in degrees)
angle = 0:binSize:180; % initialize array of angles based on bin size
minFreq = index_center * 0.1; % high pass filter: distance of pixels to NOT be counted near center of frequency image
maxFreq = index_center * 0.60; % low pass filter: distance of pixels to NOT be counted near periphery of frequency image
binSums = zeros(2,length(angle)); % first row is frequency value sum, second row is number of pixels summed over
% iterate through angles:
for a = 1:length(angle)
    % iterate through pixels in the frequency image:
    for y = index_center:length(data_freq_real-1)
        for x = 1:length(data_freq_real)
            % skip pixels if they're near the center or near the periphery (~band pass filter):
            if polarMap(y,x,1) < minFreq || polarMap(y,x,1) > maxFreq
            % check if angle of pixel is in the degree angle bin:
            elseif polarMap(y,x,2) > (angle(a) - binSize/2) && polarMap(y,x,2) <= (angle(a) + binSize/2)
                % if so, add that pixel's frequency value to the bin sum:
                binSums(1,a) = binSums(1,a) + data_freq_real(y,x);
                binSums(2,a) = binSums(2,a) + 1; % also increment number of pixels summed over
                % visualize the pixels that are being summed by blanking them out:
                data_freq_real(y,x) = 0;
            end
        end
    end
end
% normalize each bin sum by number of pixels summed over in that bin (though each bin should have ~ the same number of pixels):
binSums = binSums(1,:) ./ binSums(2,:);
binSums(isnan(binSums)) = 0; % to account for last bin, if no pixels
% baseline and normalize bin sums so that they together sum to 1 for each subdivision:
binSums = binSums - min(binSums);
%binSums = binSums / sum(binSums);
% flip order (because computed "upside-down"):
binSums = flip(binSums);
        
% visualize the pixels that were summed:
subplot(2,3,4);imshow(log(data_freq_real),[]); set(gca,'XTick',[],'YTick',[]); title('pixels counted'); % show image in subplot
% display polar plot:
subplot(2,3,5); polarplot(0:deg2rad(binSize):2*pi,[binSums,binSums(2:end)]); title('alignment frequency'); %rlim([0,0.015]); 

% orient binSums aroung max value:
[~,index_max] = max(binSums);
angle_max = angle(index_max);
binSums_aroundMax = [];
for i = -(length(binSums)-1)/2:(length(binSums)-1)/2
    if i+index_max >= 1 && i+index_max <= length(binSums)
        binSums_aroundMax(end+1) = binSums(i+index_max);
    elseif i+index_max < 1
        binSums_aroundMax(end+1) = binSums(end + (i+index_max));
    elseif i+index_max > length(binSums)
        binSums_aroundMax(end+1) = binSums(i+index_max - length(binSums));     
    end
end

% convert binSums at each angle into vectors and then average them to get mean orientation:
vector = zeros(length(angle),2); % initialize array of vectors
vector_mean = zeros(1,2);
for a = 1:length(angle)
    % calculate and store the x component of the vector:
    vector(a,1) = binSums_aroundMax(a) * cos(deg2rad(angle(a)+angle_max-90));
    % calculate and store the y component of the vector:
    vector(a,2) = binSums_aroundMax(a) * sin(deg2rad(angle(a)+angle_max-90));
end
% mean vector:
vector_mean(1) = mean(vector(:,1));
vector_mean(2) = mean(vector(:,2));
vector_mean = vector_mean / norm(vector_mean);
% plot:
subplot(2,3,6); plot(vector(:,1),vector(:,2)); hold on
plot([-vector_mean(1), vector_mean(1)],[-vector_mean(2), vector_mean(2)]);
xlim([-0.08,0.08]); ylim([-0.08,0.08])

close


% subplot(2,3,6); imshow(img_spat,[]); set(gca,'XTick',[],'YTick',[]); hold on
% vectorScale = 20;
% x_center = length(img_spat)/2;
% y_center = length(img_spat)/2;
% line([x_center + vector_mean(1)*vectorScale, x_center - vector_mean(1)*vectorScale],...
%     [y_center + vector_mean(2)*vectorScale, y_center - vector_mean(2)*vectorScale],'Color','r','LineWidth',2);

end









