function avg_activity = avg_data_curve_fitting_mixtures_v2(neural_data, time_stamps, window2use, group2use)
    groups = string(fieldnames(neural_data));
%     group2use = 1;
    avg_activity = zeros(size(neural_data.(groups(group2use))));
    for i=1:size(neural_data.(groups(group2use)),2) %loop through neurons
        for j=1:size(neural_data.(groups(group2use)),1) %loop taste stimuli
            avg_data_tmp = zeros(size(neural_data.(groups(group2use)){j,i},1),1); %avg activity for all trials of a given taste and neuron
            for p=1:size(neural_data.(groups(group2use)){j,i},1) %loop through trials
                ts2use = calculate_ts2use(time_stamps.(groups(group2use)){j,i}(p,:), window2use);
                avg_data_tmp(p) = mean(neural_data.(groups(group2use)){j,i}(p,ts2use));
            end
            avg_activity(j,i) = mean(avg_data_tmp);
        end
    end
end

