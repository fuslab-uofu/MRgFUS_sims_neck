% This function plots a slice of a pressure or Q volume, overlaid on a grayscale slice of model
function [figHandle] = plotPQSlice_neck(modelSlice, PQSliceMax, titleMax, PQSliceTarget, titleTarget, displayLimits, voxXY, figHandle)
    
    if nargin < 8
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
    imagesc(baseAxes_max, modelSlice);
    colormap(baseAxes_max, 'gray');
    axis(baseAxes_max, 'image');
    axis(baseAxes_max, 'xy');
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
    axis(overlayAxes_max, 'image');
    axis(overlayAxes_max, 'xy');
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
    
   
    % On the right, plot P/Q slice at target ----------
    % model in grayscale
    baseAxes_target = axes('Position', [0.5 0.05 0.4 0.9]);
    imagesc(baseAxes_target, modelSlice);
    colormap(baseAxes_target, 'gray');
    axis(baseAxes_target, 'image');
    axis(baseAxes_target, 'xy');
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
    axis(overlayAxes_target, 'image');
    axis(overlayAxes_target, 'xy');
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


    %%
    % plot voxel marker(s)
    hold on;
    plot(voxXY(1), voxXY(2), 'rs', 'linewidth', 2);
    axis xy;

    function syncAxesPosition(ax1, ax2)
        ax2.Position = ax1.Position;
    end
end