function save_data_decoding_mixtures_v3(responsive_neurons, frames,  neural_data, window2use, num_neurons2use, name, path2use)
    binned_data = []; binned_labels = []; binned_site_info = [];

    % neural_data2use = neural_data(neurons2use2);

    neurons2use = responsive_neurons;
%     neurons2use3 = 1:size(neural_data,2);
%     neurons2use = sort(neurons2use3(randperm(length(neurons2use3),num_neurons2use))); %randomly subsample neurnos
    qqq=1;
    for i=1:length(neurons2use)
        data2use = []; trialID_tmp = []; trialID_side_tmp=[];
        for j=1:size(neural_data, 1)
            data2use = [data2use; neural_data{j,neurons2use(i)}];
            trialID_tmp = [trialID_tmp; repmat(j,size(neural_data{j,neurons2use(i)},1),1)];
            if j<4
                r=1;
            else
                r=2;
            end
            trialID_side_tmp = [trialID_side_tmp;  repmat(r,size(neural_data{j,neurons2use(i)},1),1)];
            num_trials_taste(j,i) = size(neural_data{j,neurons2use(i)},1);
            
        end
        %use timestamps
        time_stamps = mean(frames{1,neurons2use(i)});
        %calculate index for time stamps
        ts2use = intersect(find(time_stamps>window2use(1)), find(time_stamps<window2use(2)));
        
        %format data for ndt decoding toolbox
%         if isempty(find(num_trials_taste(:,i)<10)) %omit neurons with less than 10 trials per taste
%             binned_data{1,qqq} = mean(data2use(:,ts2use),2);
%             binned_labels.trialID{1,qqq} = trialID_tmp;
%             binned_labels.direction{1,qqq} = trialID_side_tmp;
%             qqq=qqq+1;
%         end
        %use all neurons, regardless of trial count
        binned_data{1,i} = mean(data2use(:,ts2use),2);
        binned_labels.trialID{1,i} = trialID_tmp;
        binned_labels.direction{1,i} = trialID_side_tmp;
    end
    cd(path2use)
    save(append(name,'_all_neurons_decoding.mat'),'binned_data','binned_labels','binned_site_info')

end

