function [fitVals] = T4_toeLinearFitting(allData)

% "allData" is a 2D array of pooled stress-stretch results, with stretch
% for sample 1 in the 1st column, stress for sample 1 in the 2nd column,
% stretch for sample 2 in the 3rd column, etc.
% "fitVals" wil be a data structure containing E_toe, E_lin, and λ*

fitVals = zeros(size(allData,2)/2,3);

for dataSetIndex = 1:2:size(allData,2)

    % extract current data set from all data:
    data = allData(:,dataSetIndex:dataSetIndex+1);
    data = data(1:find(data(:,1),1,'last'),:); % removes trailing zeros
       
    % smooth and normalize data:
    xvals = 1:0.1:data(end,1);
    smoothingFit = fit(data(:,1),data(:,2),'power1');
    smoothData = [xvals', smoothingFit(xvals')];
    zeroedData = [smoothData(:,1)-smoothData(1,1), smoothData(:,2)-smoothData(1,2)];
    normalizedData = [zeroedData(:,1)/max(zeroedData(:,1)), zeroedData(:,2)/max(zeroedData(:,2))];
    
    % find index of elbow (defined as the point of max distance to diagonal):
    distToDiag = normalizedData(:,1) - normalizedData(:,2);
    [~,elbowIndex] = max(distToDiag);
    
    % define -/+15% x-axis bounds around elbow:
    xRange = (max(smoothData(:,1)) - min(smoothData(:,1))) * 0.15;
    xBounds = [smoothData(elbowIndex,1) - xRange, smoothData(elbowIndex,1) + xRange];
    
    % use bounds to split data into 'toe' and 'linear' subsets:
    data_toe = data((data(:,1) < xBounds(1)), :);
    data_lin = data((data(:,1) > xBounds(2)), :);
    
    % perform linear fits for the 'toe' and 'linear' subsets:
    fit_toe = polyfit(data_toe(:,1),data_toe(:,2),1);
    fit_lin = polyfit(data_lin(:,1),data_lin(:,2),1);
    
    % calculate toe and linear region moduli and transition strain:
    E_toe = fit_toe(1);
    E_lin = fit_lin(1);
    stretch_trans = (fit_toe(2)-fit_lin(2)) / (fit_lin(1)-fit_toe(1)); % horizontal value (strain) of intersect
    fitVals(1+((dataSetIndex-1)/2),:) = [E_toe,E_lin,stretch_trans];
    
    % plot data:
    figure('Position', [20+(5*dataSetIndex) 50+(5*dataSetIndex) 250 250]);
    plot(data(:,1),data(:,2),'Color','k','Linewidth',4); hold on
    plot(smoothData(:,1),smoothData(:,2),'Color','g','Linewidth',1);
    % visualize elbow zone:
    scatter(smoothData(elbowIndex,1),smoothData(elbowIndex,2),'filled','g')
    line([xBounds(1),xBounds(1)],[0,max(smoothData(:,2))],'Color','g','LineStyle','--');
    line([xBounds(2),xBounds(2)],[0,max(smoothData(:,2))],'Color','g','LineStyle','--');
%     plot(data_toe(:,1),data_toe(:,2),'Color','r','Linewidth',4);
%     plot(data_lin(:,1),data_lin(:,2),'Color','m','Linewidth',4);
    % visualize 'toe' and 'linear' region linear fits:
    plot(smoothData(:,1),polyval(fit_toe,smoothData(:,1)),'Color','r','Linewidth',2);
    plot(smoothData(:,1),polyval(fit_lin,smoothData(:,1)),'Color','m','Linewidth',2);
    line([stretch_trans,stretch_trans],[0,max(smoothData(:,2))],'Color','b');
    % make look nicer:
    ylim([0,max(smoothData(:,2))])
    xlim([1,max(smoothData(:,1))])
    % display calculated values on plot:
    text(1+max(smoothData(:,1))*0.03,max(smoothData(:,2))*0.91,['E_t = ' num2str(round(E_toe,1)) ' kPa'],'Color','r','FontSize',12);
    text(1+max(smoothData(:,1))*0.03,max(smoothData(:,2))*0.81,['E_l = ' num2str(round(E_lin,1)) ' kPa'],'Color','m','FontSize',12);
    text(1+max(smoothData(:,1))*0.03,max(smoothData(:,2))*0.74,['λ* = ' num2str(round(stretch_trans,1))],'Color','b','FontSize',12);

end

% reformat output array with column titles:
fitVals = [{'E_toe (kPa)','E_lin (kPa)','λ*'}; num2cell(fitVals)]; 

end

