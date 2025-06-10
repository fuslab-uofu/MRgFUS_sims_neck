function [figHandle] = plotTSlice(modelSlice, TSlice, myTitle, displayLimits, targetXY)
    % This function plots a slice of a temperatures volume, overlaid on a grayscale slice of model
    
    figHandle = figure; 

    % pick colormap based on which .m file exists:
    if exist('magmaColorMap.m','file')
        myColorMap = magmaColorMap();
    elseif exist('fireColorMap.m','file')
        myColorMap = fireColorMap();
    else
        % fallback if neither exists
        myColorMap = hot(256);
    end

    % plot the grayscale model
    baseAxes = axes(figHandle);
    imagesc(baseAxes, modelSlice);
    colormap(baseAxes,'gray');
    axis image;
    axis xy;

    % plot a legend with media descriptions
    legends = {'Water', 'Bone/Spine', 'Skin', 'Spinal Cord', 'Blood Vessels', 'Fat', 'Muscle', 'CSF'};
    hold on;
    for K = 1:8
        hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor','none'); %#ok<AGROW>
    end
    uistack(hidden_h,'bottom');
    legend(hidden_h, legends(1:8), 'Location', 'northwestoutside');
    tt = title(myTitle);
    tt.FontSize = 8;
    
    ylabel('Row');
    xlabel('Page');

    % the target voxel
    scatter(targetXY(1), targetXY(2), 'gs', 'linewidth', 2, 'MarkerEdgeAlpha', 0.5);

    % plot temperature slice
    overlayAxes = axes(figHandle);
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
    overlayAxes.Position = baseAxes.Position;

    cb = colorbar(overlayAxes);
    cb.Label.String = 'Temperature (°C)';
    clim(overlayAxes, [displayLimits(1), displayLimits(2)]);

    baseAxes.Position = overlayAxes.Position;
    set(overlayAxes, 'Visible', 'off');
    axis image;
    axis xy;
   
end