function [selective_neurons_same_stimulus, selective_neurons_same_choice, ...
    selective_neurons_same_stimulus_by_direction, selective_neurons_same_choice_by_direction] = ...
    find_selective_neurons_correct_error_v2(neural_data, time_stamps, trialID, resp_window)
    %error selectivity
    groups = string(fieldnames(neural_data));

    for h = 1:length(groups)
        pval_resp_same_stimulus=[];
        for z=1:size(neural_data.(groups(h)),2)
            avg_activity=[];

            for j=1:size(neural_data.(groups(h)),1)
                q=1; r=1;
                for i=1:size(neural_data.(groups(h)){j,z},1)
                    response_ts = calculate_ts2use(time_stamps.(groups(h)){j,z}(i,:),resp_window(h,:));

                    if trialID.(groups(h)){j,z}(i) == 1
                        avg_activity{1,j}(q) = mean(neural_data.(groups(h)){j,z}(i,response_ts));
                        q=q+1;
                    else
                        avg_activity{2,j}(r) = mean(neural_data.(groups(h)){j,z}(i,response_ts));
                        r=r+1;
                    end
                end

            end
            pval_resp_same_stimulus(1,z) = ranksum(avg_activity{1,1}, avg_activity{2,1}); %changed to 2
            pval_resp_same_stimulus(2,z) = ranksum(avg_activity{1,2}, avg_activity{2,2}); %changed to 2

            pval_resp_same_choice(1,z) = ranksum(avg_activity{1,1}, avg_activity{2,2}); %changed to 2
            pval_resp_same_choice(2,z) = ranksum(avg_activity{1,2}, avg_activity{2,1}); %changed to 2
            %             pval_resp(j,z) = 1;
        end
        alpha2use = 0.05/(size(neural_data.(groups(h)),1));
        %         selective_neurons.(groups(h)){1} = intersect(find(pval_resp>0),find(abs(pval_resp)<alpha2use)); %taste 1 selective
        % %         selective_neurons.(groups(h)){2} = intersect(find(pval_resp<0),find(abs(pval_resp)<alpha2use));% taste 2 selective
        %         selective_neurons.(groups(h)){1} = find(abs(pval_resp)<alpha2use); %all selective
        selective_neurons_same_stimulus.(groups(h)) = unique([find(pval_resp_same_stimulus(1,:)<alpha2use) find(pval_resp_same_stimulus(2,:)<alpha2use)]);
        selective_neurons_same_choice.(groups(h)) = unique([find(pval_resp_same_choice(1,:)<alpha2use) find(pval_resp_same_choice(2,:)<alpha2use)]);
        
        selective_neurons_same_stimulus_by_direction.(groups(h)){1} = find(pval_resp_same_stimulus(1,:)<alpha2use);
        selective_neurons_same_stimulus_by_direction.(groups(h)){2} = find(pval_resp_same_stimulus(2,:)<alpha2use);
        
        selective_neurons_same_choice_by_direction.(groups(h)){1} = find(pval_resp_same_choice(1,:)<alpha2use);
        selective_neurons_same_choice_by_direction.(groups(h)){2} = find(pval_resp_same_choice(2,:)<alpha2use);
        
%         pvals_same_choice.(groups(h)) = pval_resp_same_choice;
%         pvals_same_choice.(groups(h)) = pval_resp_same_choice;


    end
end