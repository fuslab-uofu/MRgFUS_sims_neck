% This function plots a slice of a pressure or Q volume, overlaid on a grayscale slice of model
function [figHandle] = plotPQSlice(modelSlice, PQSliceMax, titleMax, PQSliceTarget, titleTarget, displayLimits, voxXY)
   
    figHandle = figure;
    layout = tiledlayout(1, 2);

    % plot pressure at target on the left
    baseAxes_target = nexttile(layout);
    % grayscale model
    imagesc(baseAxes_target, modelSlice);
    colormap(baseAxes,'gray');
    axis image;
    axis xy;
    % plot a legend with media descriptions
    legends = {'Water', 'Bone/Spine', 'Skin', 'Spinal Cord', 'Blood Vessels', 'Fat', 'Muscle', 'CSF'};
    hold on;
    for K = 1:8
        hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor','none'); %#ok<AGROW>
    end
    hold off;
    uistack(hidden_h,'bottom');
    legend(hidden_h, legends(1:8),'Location','northwest');
    title(myTitle,'FontSize',10);
    ylabel('Row');
    xlabel('Page');
    %%
    % plot P or Q
    overlayAxes = axes('Units','normalized','position',axesPosition);
    overlayPlot = imagesc(PQSlice);
    % mask/threshold P/Q for display
    mask = ones(size(PQSlice));
    mask(PQSlice<displayLimits(1)) = 0;
    mask(PQSlice>=displayLimits(1)) = 0.7;
    set(overlayPlot,'AlphaData', squeeze(mask));
    set(overlayAxes,'Visible','off');
    colorbar(overlayAxes);
    clim(overlayAxes,[displayLimits(1),displayLimits(2)]);
    axis image;
    set(baseAxes,'Position',get(overlayAxes,'Position'));
    % plot voxel marker
    hold on;
    plot(voxXY(1), voxXY(2), 'rs', 'linewidth', 2);
    axis xy;
end