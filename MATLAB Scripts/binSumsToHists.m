function [hists_all, hist_mean] = binSumsToHists(binSums_all, binNum, centerHists)

% ensure that binNum is odd if centerHist is true:
if centerHists ~= 0 && centerHists ~= 1
    disp("ERROR: 'centerHists' must be 0 (false) or 1 (true)")
    return
end
if centerHists == 1 && rem(binNum,2) == 0
    disp("ERROR: to center histograms, 'binNum' but be an odd number")
    return
end

% initialize data variables:
binData = zeros(size(binSums_all,1), binNum);
% define bin ranges:
binSpacing = 180 / binNum;
binRanges = -90:binSpacing:90;
xVals = (binRanges(1:end-1) + binRanges(2:end)) ./ 2;
% iterate through images:
for i = 1:size(binSums_all,1)
    % convert binSums to angles for making hisograms:
    angles = [];
    for j = 1:size(binSums_all,2)
        for k = 1:round(binSums_all(i,j))
            angles(end+1) = size(binSums_all,2)-j-size(binSums_all,2)/2;
        end
    end
    % split data:
    for j = 1:binNum
        % extract data within bin ranges:
        currentBin = angles(angles>=binRanges(j));
        if j < binNum
            currentBin = currentBin(currentBin<binRanges(j+1));
        end
        % calculate summary statistics:
        binData(i,j) = length(currentBin);
    end
    % normalize so area under curve = 1:
    areaUnderCurve = sum(binData(i,:) .* binSpacing);
    binData(i,:) = binData(i,:) ./ areaUnderCurve;
    % shift to center largest bin of histogram on 0 degrees:
    if centerHists == 1
        i_zero = (binNum - 1) / 2 + 1;
        [~,i_largestBin] = max(binData(i,:));
        binData(i,:) = circshift(binData(i,:), i_zero - i_largestBin);
    end
end

% construct output:
hists_all = [xVals', binData'];
hist_mean = [xVals', mean(binData)', std(binData)'];

% show outputs as bar chart:
figure
for i = 2:size(binSums_all,1)+1
    color_bar = [0.5, 0.5, 0.5];
    bar(hists_all(:,1), hists_all(:,i), 1,'FaceColor', color_bar, 'FaceAlpha', 0.2); hold on
end
ylim([0,0.025])
xlim([-90,90])
xlabel('angle (degrees)')
ylabel('relative frequency')

figure
bar(hist_mean(:,1), hist_mean(:,2), 1,'FaceColor', [0.85, 0.85, 0.85]); hold on
ylim([0,0.025])
xlim([-90,90])
xlabel('angle (degrees)')
ylabel('relative frequency')
errorbar(hist_mean(:,1), hist_mean(:,2), hist_mean(:,3), hist_mean(:,3), 'Color', [0 0 0], 'LineStyle', 'none', 'CapSize', 2);

end





