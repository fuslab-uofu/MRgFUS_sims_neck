function [figHandle] = plotTSlice_neck(modelSlice, TSlice, myTitle, displayLimits, targetXY, otherXY, figHandle)
    % This function plots a slice of a temperatures volume, overlaid on a grayscale slice of model
    
     % Create figure at 80% of screen size
        % Get screen size in normalized units
    screenSize = get(0, 'ScreenSize');  % [left bottom width height] in pixels
    
    % Define desired figure size as a percentage of screen size
    widthPercent = 0.6;
    heightPercent = 0.6;
    
    % Calculate figure position
    figWidth = screenSize(3) * widthPercent;
    figHeight = screenSize(4) * heightPercent;
    figLeft = (screenSize(3) - figWidth) / 2;
    figBottom = (screenSize(4) - figHeight) / 2;
    figHandle = figure('Units', 'pixels', ...
                       'Position', [figLeft, figBottom, figWidth, figHeight]);
    figHandle.Units = 'normalized';

    % pick colormap based on which .m file exists:
    if exist('colorMap_magma.m','file')
        myColorMap = colorMap_magma();
    elseif exist('colorMap_fire.m','file')
        myColorMap = colorMap_fire();
    else
        % fallback if neither exists
        myColorMap = hot(256);
    end

    % plot the grayscale model
    baseAxes = axes(figHandle);
    imagesc(baseAxes, modelSlice);
    colormap(baseAxes,'gray');
    
    % plot a legend with media descriptions
    legends = {'Water', 'Bone/Spine', 'Skin', 'Spinal Cord', 'Blood Vessels', 'Fat', 'Muscle', 'CSF'};
    hold on;
    for K = 1:8
        hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor','none'); %#ok<AGROW>
    end
    uistack(hidden_h,'bottom');
    legend(hidden_h, legends(1:8), 'Location', 'northwest');
    tt = title(myTitle);
    tt.FontSize = 8;
    
    ylabel('Row');
    xlabel('Page');

    % plot temperature slice
    overlayAxes = axes('Position', baseAxes.Position);
    overlayPlot = imagesc(TSlice);
    colormap(overlayAxes, myColorMap);

    % mask/threshold for display
    mask = ones(size(TSlice));
    if (length(displayLimits) == 2)
        mask(TSlice < displayLimits(1)) = 0;
        mask(TSlice >= displayLimits(1)) = 0.9;
    elseif (length(displayLimits) > 2) % plot target volume darker
        mask(TSlice < displayLimits(1)) = 0;
        mask(TSlice >= displayLimits(1)) = 0.8; % everything
    
        ind1 = displayLimits(4):displayLimits(4 + displayLimits(3)-1);
        ind2 = displayLimits(4 + displayLimits(3)):displayLimits(end);
        mask(ind1, ind2) = 0.9; %1; % target volume
    end
    overlayPlot.AlphaData = squeeze(mask);

    cb = colorbar(overlayAxes);
    cb.Label.String = 'Temperature (°C)';
    clim(overlayAxes, [displayLimits(1), displayLimits(2)]);

    baseAxes.Position = overlayAxes.Position;
    set(overlayAxes, 'Visible', 'off');

    % plot target voxel point
    hold on;
    plot(targetXY(1, 1), targetXY(1, 2), 'r', 'Marker', '*'); 
    % plot other points of interest
    scatter(otherXY(:, 1), otherXY(:, 2), 25, 'k');

    axis(baseAxes, 'image');
    axis(baseAxes, 'xy');
    axis(overlayAxes, 'image');
    axis(overlayAxes, 'xy');
   
end