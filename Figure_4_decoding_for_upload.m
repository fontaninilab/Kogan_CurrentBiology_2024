%% decoding analysis using neural decoding tool box on taste mixtures in sampling and delay epochs
%this is only looking at correct trials, same as previous analysis
%
%groups for analysis, central and lateral alignments
%central used for sampling epoch and lateral used for delay
%central aligned to first lick to central spout (sampling onset)
%lateral aligned to first lick to lateral (choice onset)
%clear all
clear; clc; close all
%add current folder and subfolders containing utility functions and data
%files
addpath(genpath(pwd))
load('neural_data_mixtures.mat')
path2use = fullfile(pwd,'decoding_data');
if ~isfolder(path2use)
    mkdir(path2use)
end
%%
groups = string(fieldnames(neural_data)); 

%parameters for decoding
num_repeats = 1; 
num_neurons2use = 2000; %how many neurons for pseudopopulation
%use all tastes
tastes2use = [1 2 3 4 5 6];
shuffle = 0;

%use mixture pairs
tastes2use_pairs = [1 6;2 5; 3 4] ;

%response windows used, can also run it on the baseline windows and
%performance drops 
resp_window_groups = [0 1.5; -1.5 0; 0 1.5; -1.5 0];
baseline_window = [-4 -2.5; -6.5 -5; -4 -2.5; -6.5 -5];

confusion=zeros(length(tastes2use),length(tastes2use),length(groups)) ;
performance=zeros(1,length(groups));
performance_pairs = zeros(size(tastes2use_pairs,1), length(groups));

%run decoding on all mixtures
for h=1:length(groups) %loop through 4 groups, sampling/delay for pre and post
    group2use = h;
    
    %save the data, this is ignoring the responsive neurons input, but can
    %be changed to use only those neurons, it will instead randomly
    %subsample neurons from the total population using the num_neurons2use
    %variable, 
    %
    %data will be saved in path2use folder
    
    %use same group of neurons for both sampling and delay
    if ismember(h,[1 3])
        neurons2use=[]; all_neurons=[];
        all_neurons = 1:size(neural_data.(groups(group2use)),2);
        neurons2use = sort(all_neurons(randperm(length(all_neurons),num_neurons2use))); %randomly subsample neurons
    end
    save_data_decoding_mixtures_v3(neurons2use, time_stamps.(groups(group2use)),...
        neural_data.(groups(group2use)),resp_window_groups(group2use, :), num_neurons2use, char((groups(group2use))),path2use)
    
    %which data to load
    binned_format_file_name = [char((groups(group2use))),'_all_neurons_decoding.mat'];
    
    %decoding performance on all mixtures
    [performance(h),decoding] = ndt_decoding_mixtures_v2(binned_format_file_name,shuffle, tastes2use,'trialID', num_neurons2use);
    confusion(:,:,h) = decoding.NORMALIZED_RANK_RESULTS.confusion_matrix_results.rank_confusion_matrix;
    
    %decoding on mixture pairs using the same neurons as above
    for p=1:size(tastes2use_pairs,1)
        [performance_pairs(p,h),~] = ndt_decoding_mixtures_v2(binned_format_file_name,shuffle, tastes2use_pairs(p,:),'trialID', num_neurons2use);
        
    end
end

%% plot Figure 4
lbls = ["100","75","60","40","25","0"];
titles = ["Pre Sampling","Pre Delay","Post Sampling","Post Delay"];
figure
%plot performance for mixture pairs
subplot(2,2,1)
hold on
plot([3 2 1],mean(performance_pairs(:,1),2),'linewidth',3,'color','k','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor','k')
plot([3 2 1],mean(performance_pairs(:,2),2),'linewidth',3,'color','c','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor','c')

xlim2use = get(gca, 'xlim');
line(xlim2use, [.5 .5], 'linestyle','--','color','k')
ylim([.4 1])
% yticks(.4:.1:1)
xticks([1 2 3])
xticklabels(["60 vs 40", "75 vs 25", "100 vs 0"]) 
xlabel('Sucrose (mM)')
ylabel('Performance')
title('Decoding - Mixture Pairs')
box off
axis square
legend('Sampling','Delay','Location','east')

%

%plot overall performance for all mixtures
subplot(2,2,2)
b1 = bar([1.3 1.7], performance([1 2]),'BarWidth',.6);
hold on
b1.FaceColor = 'flat';
line([1 2],[1/6 1/6],'linestyle','--','color','k')
b1.CData(1,:) = [0 0 0]; 
b1.CData(2,:) = [0 1 1]; 
xticks([1.3 1.7])
xlim([1 2])
xticklabels(["Sampling","Delay"])
yticks([0 .3 .6])
ylim([0 .6])
axis square
box off
ylabel('Performance')
title('Mixture Decoding Performance')

%plot confusion
for h=1:2
    subplot(2,2,h+2)
    
    imagesc(1:6, 1:6, confusion(:,:,h), [.3 .85])
    cb = colorbar;
    cb.Ticks = [.3 .85];
    psn_tmp = cb.Position;
    cb.Position = (psn_tmp.*[1 1 1 .6])+[.06 0 0 0];
    ylabel(cb,'Performance')
    cb.Label.Position = cb.Label.Position - [1.5 0 0];
    box off
    axis square
    xticks(1:6)
    yticks(1:6)
    xticklabels(lbls)
    yticklabels(lbls)
    title(titles(h))
    xlabel('Sucrose (mM)')
    ylabel('Sucrose (mM)')

end
