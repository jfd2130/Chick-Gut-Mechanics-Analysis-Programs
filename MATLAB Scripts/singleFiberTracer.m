function [rawTrace, smoothTrace, straightness, direction] = singleFiberTracer()
%% Setup Image

% I. allow the user to select an image file:
[fileName, pathName] = uigetfile('*.tif', 'Open SHG image (tif):');
image = imread(fullfile(pathName, fileName));

% II. display the buttons and image:
fig_finishButton = figure(1);
buttonHandle_finishTrace = uicontrol('Style', 'PushButton','String', 'Finish Trace','BackgroundColor','green','FontWeight','bold','Callback','delete(gcbf)');
buttonHandle_finishTrace.Position = [10 10 100 60]; % [x pos frombottom left, y pos from bottom left, width, height]
fig_finishButton.Position = [300 700 120 80];

fig_resetButton = figure(2);
buttonHandle_resetTrace = uicontrol('Style', 'PushButton','String', 'Reset Trace','BackgroundColor','yellow','FontWeight','bold','Callback','delete(gcbf)');
buttonHandle_resetTrace.Position = [10 10 100 60]; % [x pos frombottom left, y pos from bottom left, width, height]
fig_resetButton.Position = [300 500 120 80];

fig_endButton = figure(3);
buttonHandle_endTracing = uicontrol('Style', 'PushButton','String', 'End Tracing','BackgroundColor','red','FontWeight','bold','Callback','delete(gcbf)');
buttonHandle_endTracing.Position = [10 10 100 60];
fig_endButton.Position = [300 300 120 80];

figure(4);
imshow(image); hold on
set(gcf, 'units', 'normalized', 'position', [0.25 0 0.75 0.925])

% III. make variables to hold outputs:
numFibers_max = 1000; % max number of fibers that can be traced from one image
numNodes_max = 100; % max number of nodes that can be in one fiber
numSegsInBezier = 11; % how many segments are calculated between each node when smoothing with a Bezier curve
rawTrace = zeros(numFibers_max, numNodes_max, 2);
smoothTrace = zeros(numFibers_max, numNodes_max * numSegsInBezier, 2);
straightness = zeros(numFibers_max, 1); % one calculated per fiber (using distances between points in smooth trace)
direction = zeros(numFibers_max, 1); % one calculated per fiber (using linear regression through points in smooth trace)


%% Manual Fiber Tracing:

% tracing iterations:
for fiberID = 1:numFibers_max
    tracingFinished = false; % used for cutting the while loop
    while ~tracingFinished
        % initialize a variable to hold displayed points:
        scatterPoints = zeros(numNodes_max, 1);
        % node selection:
        for nodeID = 1:numNodes_max
            % collect one node position:
            [nodeX, nodeY] = ginputWhite(1); % ginputWhite() is a function you should be able get online, otherwise just replace with ginput(()
            % check if any button has already been pressed:
            if ~ishandle(buttonHandle_finishTrace) && nodeID > 1
                % reset button figure:
                fig_finishButton = figure(1);
                buttonHandle_finishTrace = uicontrol('Style', 'PushButton','String', 'Finish Trace','BackgroundColor','green','FontWeight','bold','Callback','delete(gcbf)');
                buttonHandle_finishTrace.Position = [10 10 100 60]; % [x pos frombottom left, y pos from bottom left, width, height]
                fig_finishButton.Position = [300 700 120 80];
                set(gca,'visible','off')
                % calculate Bezier curve to interpolate points between nodes: (SOURCE: Aaron Wetzler, aaronwetzler@gmail.com, 2009)
                [controlPoint1, controlPoint2] = findControlPoints(squeeze(rawTrace(fiberID, 1:nodeID-1, :)));
                bezierCurve = getBezier(squeeze(rawTrace(fiberID, 1:nodeID-1, :)), controlPoint1, controlPoint2);
                smoothTrace(fiberID, 1:length(bezierCurve), :) = bezierCurve; % store 'smoothed' points
                figure(4) % bring image figure back to focus
                delete(scatterPoints); % get rid of the displayed selected points
                plot(bezierCurve(:, 1), bezierCurve(:, 2), 'r'); % display the curve
                % end tracing and inner for loop:
                disp(['fiber ' num2str(fiberID) ' finished'])
                tracingFinished = true; % will stop while loop
                break;
            elseif ~ishandle(buttonHandle_resetTrace)
                % reset button figure:
                fig_resetButton = figure(2);
                buttonHandle_resetTrace = uicontrol('Style', 'PushButton','String', 'Reset Trace','BackgroundColor','yellow','FontWeight','bold','Callback','delete(gcbf)');
                buttonHandle_resetTrace.Position = [10 10 100 60]; % [x pos frombottom left, y pos from bottom left, width, height]
                fig_resetButton.Position = [300 500 120 80];
                set(gca,'visible','off')
                % reset trace for current fiber:
                rawTrace(fiberID, :, :) = 0; % clear the node positions for this fiber:
                figure(4) % bring image figure back to focus
                delete(scatterPoints); % get rid of the displayed selected points
                % end inner for loop only (allow this fiber to be retraced):
                disp('tracing reset')
                break;
            elseif ~ishandle(buttonHandle_endTracing)
                close all
                % end tracing and outer for loop:
                disp('Fiber tracing ended. Remember to save your data!');
                tracingFinished = true; % will stop while loop
                break;
            else
                % if no buttons were pressed, store the node position to the current trace:
                rawTrace(fiberID, nodeID, 1) = nodeX;
                rawTrace(fiberID, nodeID, 2) = nodeY;
                % draw new node onto image:
                scatterPoints(nodeID) = scatter(nodeX, nodeY, 150, 'g', '.');
            end
        end
    end
    if ~ishandle(buttonHandle_endTracing)
        % end outer for loop:
    	break;
    end
end


%% Calculate straightness and direction:

for fiberID = 1:numFibers_max
    % do not proceed if we've reached the last traced fiber (based on this fiber's first node's (x) position being 0):
    if smoothTrace(fiberID, 1, 1) == 0
        break;
    end
    
    % calculate straightness as end-to-end length / contour length:
    length_ends = 0;
    length_contour = 0;
    for nodeID = 2:(numNodes_max * numSegsInBezier)
        % do not proceed if we've reached the last node (based on this node's (x) position being 0):
        if smoothTrace(fiberID, nodeID, 1) == 0
            nodeID_last = nodeID - 1; % store for convinent use later
            % use pathagorian theorum to calculate end-to-end length:
            length_ends = sqrt(sum((smoothTrace(fiberID, nodeID - 1, :) - smoothTrace(fiberID, 1, :)).^2));
            % exit loop:
            break;
        end
        % use pathagorian theorum to add distance between current and previous node:
        length_contour = length_contour + sqrt(sum((smoothTrace(fiberID, nodeID, :) - smoothTrace(fiberID, nodeID - 1, :)).^2));
    end
    straightness(fiberID) = length_ends / length_contour; % store straightness value
    
    % calculate direction as the angle in radians from 0 (in [pi/2,-pi/2]) of a line of regression:
    regression = polyfit(smoothTrace(fiberID, 1:nodeID_last, 1), smoothTrace(fiberID, 1:nodeID_last, 2), 1);  % linear regression, regression(1) gives slope
    direction(fiberID) = -atan2(regression(1), 1); % negative b/c image y coordinates are "upside down"
end
pause(5)
end








