function [binData] = bin2Ddata(data, binAxis, binNum)
% 'data' should be Nx2, where N is total # data points
% 'binAxis' should be 1 or 2
% 'binNum' should be an integer > 0

if (binAxis ~= 1 && binAxis ~= 2)
    disp('Unknown axis: Please enter only "1" or "2" for binAxis.');
    return 
end
binNum = round(abs(binNum)); % in case the input was not an integer
binData = zeros(binNum, 4); % for each bin: axis 1 mean, axis 2 mean, axis 1 std, axis 2 std

% define bin ranges:
binMax = max(data(:,binAxis));
binMin = min(data(:,binAxis));
binSpacing = abs(binMax - binMin) / binNum;
binRanges = binMin:binSpacing:binMax;
% split data:
for i = 1:binNum
    % extract data within bin ranges:
    currentBin = data(data(:,binAxis)>=binRanges(i),:);
    currentBin = currentBin(currentBin(:,binAxis)<binRanges(i+1),:);
    if i == binNum
        % include the max data point in the last bin:
        currentBin = [currentBin;data(data(:,binAxis)==binMax,:)];
    end
    % calculate summary statistics:
    binData(i,1) = mean(currentBin(:,1)); % axis 1 mean
    binData(i,2) = mean(currentBin(:,2)); % axis 2 mean
    binData(i,3) = std(currentBin(:,1)); % axis 1 std
    binData(i,4) = std(currentBin(:,2)); % axis 2 std
end

end