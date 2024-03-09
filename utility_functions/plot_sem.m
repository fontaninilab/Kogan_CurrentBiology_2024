function plot_sem(x_vals, activity_matrix)
    xx = [x_vals, fliplr(x_vals)];
    yy = [mean(activity_matrix) + find_sem((activity_matrix)),...
        fliplr(mean(activity_matrix) - find_sem((activity_matrix)))];
    h = fill(xx, yy, [0.5 0.5 0.5]); % plot this first for overlay purposes
    set(h, 'facealpha', 0.3, 'edgecolor', 'none');
end