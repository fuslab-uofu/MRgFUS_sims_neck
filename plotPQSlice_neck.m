% =============================================================================
%               plotPQSlice_neck()
% -----------------------------------------------------------------------------
%
% Description:  Plots pressures or Q at the specified slice overlayed on the 
%               corresponding slice of grayscale segmented model.
%
% Authors:      Michelle Kline & Marta M. Iversen  
%               Department of Radiology and Imaging Sciences  
%               University of Utah  
%
% Inputs:       1. slice of model at max P/Q voxel
%               2. slice of P/Q at max P/Q voxel
%               3. title for max P/Q plot
%               4. slice of model at target P/Q voxel
%               5. slice of P/Q at target P/Q voxel
%               6. title for target P/Q plot
%               7. display limits [lower, upper] for colormap and masking
%               8. x y coordinate of max voxel
%               9. x y coordinate of target voxel
%               
% Output:       figure with max P/Q on the left and target P/Q on the right
%
% Requirements: MATLAB R2020a or later (toolboxes?)  
%
% Dependencies: 
%               
% =============================================================================

function [figHandle] = plotPQSlice_neck(modelSliceMax, PQSliceMax, titleMax, modelSliceTarget, PQSliceTarget, titleTarget, displayLimits, maxXY, targXY, otherXY, figHandle)
    
    if nargin < 11
        % Create figure at 80% of screen size
        % Get screen size in normalized units
        screenSize = get(0, 'ScreenSize');  % [left bottom width height] in pixels
        
        % Define desired figure size as a percentage of screen size
        widthPercent = 0.8;
        heightPercent = 0.8;
        
        % Calculate figure position
        figWidth = screenSize(3) * widthPercent;
        figHeight = screenSize(4) * heightPercent;
        figLeft = (screenSize(3) - figWidth) / 2;
        figBottom = (screenSize(4) - figHeight) / 2;
        figHandle = figure('Units', 'pixels', ...
                       'Position', [figLeft, figBottom, figWidth, figHeight]);
    end

    figHandle.Units = 'normalized';
        
    % On the left, plot P/Q slice at max ---------
    % model in grayscale
    baseAxes_max = axes('Position', [0.05 0.05 0.4 0.9]);
    imagesc(baseAxes_max, modelSliceMax);
    colormap(baseAxes_max, 'gray');
    title(baseAxes_max, titleMax);
    % plot a legend with media descriptions
    legends = {'Water', 'Bone/Spine', 'Skin', 'Spinal Cord', 'Blood Vessels', 'Fat', 'Muscle', 'CSF'};
    hold on;
    for K = 1:8
        hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor','none'); %#ok<AGROW>
    end
    uistack(hidden_h,'bottom');
    legend(hidden_h, legends(1:8),'Location','northwest');
    xlabel(baseAxes_max, 'Page');
    ylabel(baseAxes_max, 'Row');

    % overlay P or Q
    overlayAxes_max = axes('Position', baseAxes_max.Position);
    overlayPlot_max = imagesc(overlayAxes_max, PQSliceMax);
    
    overlayAxes_max.Visible = 'off';
    
    % alpha mask
    mask = ones(size(PQSliceMax));
    mask(PQSliceMax < displayLimits(1)) = 0;
    mask(PQSliceMax >= displayLimits(1)) = 0.7;
    set(overlayPlot_max, 'AlphaData', mask);
    
    % colorbar
    cb = colorbar(overlayAxes_max);
    cb.Label.String = 'Pressure (Pa)';
    clim(overlayAxes_max, [displayLimits(1), displayLimits(2)]);
    baseAxes_max.Position = overlayAxes_max.Position;
    
    % for zooming and panning
    linkaxes([baseAxes_max, overlayAxes_max]);
    % for ensuring the both base and overlay axes stay aligned when interactive resizing of figure
    figHandle.SizeChangedFcn = @(src, event) syncAxesPosition(baseAxes_max, overlayAxes_max);

    % plot max voxel point
    hold on;
    plot(maxXY(1, 1), maxXY(1, 2), 'r', 'Marker', '*'); 
    % plot other points of interest
    scatter(otherXY(:, 1), otherXY(:, 2), 25, 'k');

    axis(baseAxes_max, 'image');
    axis(baseAxes_max, 'xy');
    axis(overlayAxes_max, 'image');
    axis(overlayAxes_max, 'xy');
   
    % On the right, plot P/Q slice at target ----------
    % model in grayscale
    baseAxes_target = axes('Position', [0.5 0.05 0.4 0.9]);
    imagesc(baseAxes_target, modelSliceTarget);
    colormap(baseAxes_target, 'gray');
    title(baseAxes_target, titleTarget);
    % plot a legend with media descriptions
    legends = {'Water', 'Bone/Spine', 'Skin', 'Spinal Cord', 'Blood Vessels', 'Fat', 'Muscle', 'CSF'};
    hold on;
    for K = 1:8
        hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor','none');
    end
    uistack(hidden_h,'bottom');
    legend(hidden_h, legends(1:8),'Location','northwest');
    xlabel(baseAxes_target, 'Page');
    ylabel(baseAxes_target, 'Row');

    % overlay P or Q
    overlayAxes_target = axes('Position', baseAxes_target.Position);
    overlayPlot_target = imagesc(overlayAxes_target, PQSliceTarget);
    overlayAxes_target.Visible = 'off';
    
    % alpha mask
    mask = ones(size(PQSliceTarget));
    mask(PQSliceTarget < displayLimits(1)) = 0;
    mask(PQSliceTarget >= displayLimits(1)) = 0.7;
    set(overlayPlot_target, 'AlphaData', mask);
    
    % colorbar
    cb = colorbar(overlayAxes_target);
    cb.Label.String = 'Pressure (Pa)';
    clim(overlayAxes_target, [displayLimits(1), displayLimits(2)]);
    baseAxes_target.Position = overlayAxes_target.Position;
    
    % for zooming and panning
    linkaxes([baseAxes_target, overlayAxes_target]);
    % for ensuring the both base and overlay axes stay aligned when interactive resizing of figure
    figHandle.SizeChangedFcn = @(src, event) syncAxesPosition(baseAxes_target, overlayAxes_target);

    % plot target voxel
    hold on;
    plot(targXY(1, 1), targXY(1, 2), 'r', 'Marker', '*'); 
    scatter(otherXY(:, 1), otherXY(:, 2), 25, 'k');

    axis(baseAxes_target, 'image');
    axis(baseAxes_target, 'xy');
    axis(overlayAxes_target, 'image');
    axis(overlayAxes_target, 'xy');


    %%
    function syncAxesPosition(ax1, ax2)
        ax2.Position = ax1.Position;
    end
end