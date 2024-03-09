function [selective_neurons, pvals_all] = find_selective_neurons_correct_error(neural_data, time_stamps, trialID, resp_window)

    groups = string(fieldnames(neural_data));

    for h = 1:length(groups)
        pval_resp=[];
        for z=1:size(neural_data.(groups(h)),2)
            avg_activity=[];
            
            for j=1:size(neural_data.(groups(h)),1)
                q=1;
                for i=1:size(neural_data.(groups(h)){j,z},1)
                    if trialID.(groups(h)){j,z}(i) == 1
                    response_ts = calculate_ts2use(time_stamps.(groups(h)){j,z}(i,:),resp_window(h,:));

                    avg_activity{j}(q) = mean(neural_data.(groups(h)){j,z}(i,response_ts));
                    q=q+1;
                    end
                end

            end
            if mean(avg_activity{1}) >= mean(avg_activity{2})
                pval_resp(z) = ranksum(avg_activity{1}, avg_activity{2}); %changed to 2
            else
                pval_resp(z) = -ranksum(avg_activity{1}, avg_activity{2});
            end
        
%             pval_resp(j,z) = 1;
        end
        alpha2use = 0.05;%/(size(neural_data.(groups(h)),1));
        selective_neurons.(groups(h)){1} = intersect(find(pval_resp>0),find(abs(pval_resp)<alpha2use)); %taste 1 selective
        selective_neurons.(groups(h)){2} = intersect(find(pval_resp<0),find(abs(pval_resp)<alpha2use));% taste 2 selective
        selective_neurons.(groups(h)){3} = find(abs(pval_resp)<alpha2use); %all selective
        pvals_all.(groups(h)) = abs(pval_resp);

    end
end