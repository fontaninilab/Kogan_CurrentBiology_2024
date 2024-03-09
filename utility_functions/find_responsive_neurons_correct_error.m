function [responsive_neurons, responsive_neurons_indv_tastes, pvals_all] = ...
    find_responsive_neurons_correct_error(neural_data, time_stamps, trialID, baseline_window, resp_window)

groups = string(fieldnames(neural_data));

for h = 1:length(groups)
    for z=1:size(neural_data.(groups(h)),2)
        for j=1:size(neural_data.(groups(h)),1)
            avg_activity=[];
            q=1;
            for i=1:size(neural_data.(groups(h)){j,z},1)
                if trialID.(groups(h)){j,z}(i) == 1
                    
                baseline_ts = calculate_ts2use(time_stamps.(groups(h)){j,z}(i,:),baseline_window(h,:));
                response_ts = calculate_ts2use(time_stamps.(groups(h)){j,z}(i,:),resp_window(h,:));

                avg_activity(q,1) = mean(neural_data.(groups(h)){j,z}(i,baseline_ts));
                avg_activity(q,2) = mean(neural_data.(groups(h)){j,z}(i,response_ts));
                
                q = q+1;
                end

            end
            if diff(mean(avg_activity))>=0
                pval_resp(j,z) = ranksum(avg_activity(:,1), avg_activity(:,2));
            else
                pval_resp(j,z) = 1;
            end
        end
    end
    alpha2use = 0.05/(size(neural_data.(groups(h)),1));
%     alpha2use = 0.01;
    responsive_neurons.(groups(h)) = find(min(pval_resp)<alpha2use); %find if neuron is responsive to at least one taste
    pvals_all.(groups(h)) = pval_resp;

    for i=1:size(pval_resp,1)
        responsive_neurons_indv_tastes.(groups(h)){i} = find(pval_resp(i,:)<alpha2use);
    end
end
end