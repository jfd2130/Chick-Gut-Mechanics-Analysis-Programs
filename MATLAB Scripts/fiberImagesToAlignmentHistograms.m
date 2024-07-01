function [hists,hists_averaged] = fiberImagesToAlignmentHistograms(binNum, imgsToAvrg)
% this function applies the fiberAlignment() and binSumsToHists() functions
% to a folder of tif images (SHG max projections) and outputs each alignment
% histogram (plus an averaged UNcentered histogram):

% access image folder:
[fileName, pathName] = uigetfile('*.tif', 'Select a fiber image:');
fullName = fullfile(pathName,fileName); % these four lines are from Meshoo on MATLAB Answers...
filelist = dir([fileparts(fullName) filesep '*.tif']);
fileNames = {filelist.name}';
fileNames = fileNames(1:1:length(fileNames));
num_img = numel(filelist);
warning('off','images:initSize:adjustingMag'); % reduces Command Window clutter

% load all images into a single array:
dims_min = [1 1] .* 10e10; % used for truncating all images to the same dimensions
for img = 1:num_img
    img_current = imread(fullfile(pathName,fileNames{img}));
    if size(img_current, 1) < dims_min(1)
        dims_min(1) = size(img_current, 1);
    end
    if size(img_current, 2) < dims_min(2)
        dims_min(2) = size(img_current, 2);
    end
end
images = zeros(num_img, dims_min(1), dims_min(2));
for img = 1:num_img
    img_current = imread(fullfile(pathName,fileNames{img}));
    images(img,:,:) = img_current(1:dims_min(1), 1:dims_min(2));
end

% perform FFT fiber alignment to get bin sums for each image:
binSums = [];
for img = 1:num_img
    [binSum_current] = fiberAlignment(squeeze(images(img,:,:)));
    binSums = [binSums; binSum_current];
end

% convert bin sums to angle histograms (to center histograms, change last parameter to = 1):
[hists,~] = binSumsToHists(binSums, binNum, 0);

% average histograms for given images:
hists_averaged = hists(:,1);
for i = 1:size(imgsToAvrg, 1)
    hists_averaged(:,1+i) = zeros(size(hists_averaged,1), 1);
    num_toAvrg = 0;
    for j = 1:size(imgsToAvrg, 2)
        if imgsToAvrg(i,j) ~= 0
            hists_averaged(:,1+i) = hists_averaged(:,1+i) + hists(:,1+imgsToAvrg(i,j));
            num_toAvrg = num_toAvrg + 1;
        end
    end
    hists_averaged(:,1+i) = hists_averaged(:,1+i) ./ num_toAvrg;
end


end




