%% Scripts to plot and analyze neural data for Figure 2 of taste mixture discrimination manuscript 
% Notes: 
% 1)Load in mixture data from correct trials, pre and post, aligned to both
% central and lateral
% 2)Find responsive/selective neurons
% 3)Plot heatmaps of responsive neurons, population avg, example traces and
% population avg selectivity
% 4)Statistics and other pieces of information for the manuscript

%% Load data - correct trials
%also contains responsive/selective neurons
% cd('D:\Code_data_analysis\Curr_bio_data_upload')
%clear all
clear; clc; close all
%add current folder and subfolders containing utility functions and data
%files
addpath(genpath(pwd))
load('Figure_2_neural_data_w_significant_neurons.mat')
time_stamps = time_stamps_v2; %rename timestamps variable
%groups for organzing alignment(sampling/choice) and conditions(pre/post)
groups = string(fieldnames(neural_data_v2));
%% Calculate proportion of responsive/selective neurons
%responsive out of total
for i=1:length(groups)
   proportion_responsive_v2(i) = length(responsive_neurons_v2.(groups(i)))/size(neural_data_v2.(groups(i)),2);
   tmp(i,:) = [length(responsive_neurons_v2.(groups(i))) size(neural_data_v2.(groups(i)),2)];
end

%selective out of responsive
selectivity_numbers=[];
for i=1:length(groups)
   neurons2use = intersect(responsive_neurons_v2.(groups(i)), selective_all.(groups(i)));
   proportion_selective(i) = length(neurons2use)/size(responsive_neurons_v2.(groups(i)),2);
   selectivity_numbers(1,i) = length(neurons2use);
   selectivity_numbers(2,i) = size(responsive_neurons_v2.(groups(i)),2);
end
selectivity_numbers
%% plot example selective neurons
figure
subplot(1,2,1)
hold on
neuron2use = 1139; xlim2use = [-.5 2]; group2use = 1;
plt_example_neuron_PSTH(neuron2use, time_stamps, neural_data_v2, group2use, xlim2use)

% delay example neuron
subplot(1,2,2)
hold on
neuron2use = 583; xlim2use = [-2 .5]; group2use = 2;
plt_example_neuron_PSTH(neuron2use, time_stamps, neural_data_v2, group2use, xlim2use)

%% plot pop average all responsive neurons
figure
subplot(2,1,1)
hold on
group2use=1; xlim2use = [-1 3]; delta = 0;
neurons2use = responsive_neurons_v2.(groups(group2use));
plt_pop_PSTH(neurons2use, time_stamps, neural_data_v2, group2use, xlim2use,delta)
% ylim(ylim2use)
% yticks(0:0.05:0.15)
title(['Sampling Responsive',newline,'Pop PSTH'])

%plot pop PSTH delay
group2use=2; xlim2use = [-3 1]; delta = 0; neurons2use=[];
neurons2use = responsive_neurons_v2.(groups(group2use));
subplot(2,1,2)
hold on
plt_pop_PSTH(neurons2use, time_stamps, neural_data_v2, group2use, xlim2use,delta)
% ylim(ylim2use)
% yticks(0:0.05:0.15)
title(['Delay Responsive',newline,'Pop PSTH'])


%% plot selectivity pop average and calculate peaks for selectivity
figure
subplot(1,3,1)
hold on
group2use=1; xlim2use = [-1 3]; resp_window2use = [0 1.5]; delta = 1;
neurons2use = intersect(responsive_neurons_v2.(groups(group2use)), selective_neurons.(groups(group2use)));
[peak_delta_PSTH(group2use), peak_delta_ts(group2use)] = ...
    plt_pop_PSTH_v2(neurons2use, time_stamps, neural_data_v2, group2use, xlim2use, resp_window2use, delta);
title(['Sampling Selective',newline,'\DeltaPop PSTH'])
ylabel('\Delta\DeltaF/f')

%plot pop PSTH delay
subplot(1,3,2)
hold on
group2use=2; xlim2use = [-3 1]; resp_window2use = [-1.5 0]; delta = 1;
neurons2use = intersect(responsive_neurons_v2.(groups(group2use)), selective_neurons.(groups(group2use)));
[peak_delta_PSTH(group2use), peak_delta_ts(group2use)] = ...
    plt_pop_PSTH_v2(neurons2use, time_stamps, neural_data_v2, group2use, xlim2use, resp_window2use, delta);
title(['Delay Selective',newline,'\DeltaPop PSTH'])
ylabel('\Delta\DeltaF/f')

subplot(1,3,3)
bar([1 2],proportion_selective([1 2]),'k','BarWidth',.6)
text(1, proportion_selective(1)+.03, [num2str(selectivity_numbers(1,1)),'/',num2str(selectivity_numbers(2,1))], 'fontsize',6,'HorizontalAlignment','center')
text(2, proportion_selective(2)+.03, [num2str(selectivity_numbers(1,2)),'/',num2str(selectivity_numbers(2,2))], 'fontsize',6,'HorizontalAlignment','center')

box off
axis square
xlim([0 3])
xticklabels(["Sampling","Delay"])
ylim([0 .45])
yticks(0:.1:.4)
title(['Proportion',newline, 'Selective'])
ylabel('Proportion')

peak_delta_ts
%% Plotting heatmaps for responsive neurons

% %sampling responsive 
resp_tmp1=[]; resp_tmp2=[]; sampling_diff=[];
resp_window = [0 1.5];
whole_window = [-1 3];
group2use = 1;

neurons2use = responsive_neurons_v2.(groups(group2use));

for i=1:length(neurons2use)
    ts2use = calculate_ts2use(mean(time_stamps_v2.(groups(group2use)){1,neurons2use(i)}), resp_window);

    resp_tmp1(i,:) = mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)});
    resp_tmp2(i,:) = mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)});
    
    avg_resp_tmp1 = mean(resp_tmp1(i, ts2use));
    avg_resp_tmp2 = mean(resp_tmp2(i, ts2use));
    
    sampling_diff(i) = (avg_resp_tmp1 - avg_resp_tmp2);

end
ts2use=[];
ts2use = calculate_ts2use(mean(time_stamps_v2.(groups(group2use)){1,1}), whole_window);
resp_tmp3=[];
resp_tmp3 = [resp_tmp1(:,ts2use) resp_tmp2(:,ts2use)];

sorted_idx = [];
[~,sorted_idx] = sort(sampling_diff);

% [~,sorted_idx2] = max(resp_tmp1,[],2);
% [~,temporal_sort] = sort(sorted_idx2, 'ascend');

rowmin = min(resp_tmp3,[],2);
rowmax = max(resp_tmp3,[],2);
Brow=[];
Brow = rescale(resp_tmp3,0,1,'InputMin',rowmin,'InputMax',rowmax);

time_stamps2use = mean(time_stamps_v2.(groups(group2use)){1,1});
colorlim2use = [0 1];
% map2use = brewermap(100,'BrBG');

normalized_activity_taste1 = Brow(sorted_idx,1:(size(Brow,2)/2));
normalized_activity_taste2 = Brow(sorted_idx,1+(size(Brow,2)/2):end);

norm_diff = normalized_activity_taste1 - normalized_activity_taste2;
ts2use2 = calculate_ts2use(time_stamps2use(ts2use), resp_window);

% [~,sorted_norm_idx] = sort(mean(norm_diff(:,ts2use2),2));

%heatmap sampling responsive
figure
subplot(1,2,1)
imagesc(time_stamps2use(ts2use), 1:size(resp_tmp1,1), normalized_activity_taste1, colorlim2use)
hold on
ylim2use = [0 size(resp_tmp2,1)];
line([0 0],ylim2use, 'linestyle','--','color','k')
cb = colorbar('TickDirection','out','Ticks',[-1 0 1]);
ylabel(cb, 'Normalized \DeltaF/f');
box off
% axis square
title('Sucrose Trials')
ylabel('Neurons')
xlabel('Time (sec)')
xlim(whole_window)
yticks(round(size(resp_tmp2,1),-1))
ylim(ylim2use)

subplot(1,2,2)
imagesc(time_stamps2use(ts2use), 1:size(resp_tmp2,1),normalized_activity_taste2 , colorlim2use)
hold on
ylim2use = [0 size(resp_tmp2,1)];
line([0 0],ylim2use, 'linestyle','--','color','k')
cb = colorbar('TickDirection','out','Ticks',[-1 0 1]);
ylabel(cb, 'Normalized \DeltaF/f');
box off
% axis square
title('NaCl Trials')
ylabel('Neurons')
xlabel('Time (sec)')
xlim(whole_window)
ylim(ylim2use)
yticks(round(size(resp_tmp2,1),-1))
sgtitle('Sampling Responsive Neurons','fontsize',10)


resp_tmp1=[]; resp_tmp2=[]; sampling_diff=[];

%delay responsive
resp_window = [-1.5 0];
whole_window = [-3 1];
group2use = 2;

neurons2use = responsive_neurons_v2.(groups(group2use));

for i=1:length(neurons2use)
    ts2use = calculate_ts2use(mean(time_stamps_v2.(groups(group2use)){1,neurons2use(i)}), resp_window);

    resp_tmp1(i,:) = mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)});
    resp_tmp2(i,:) = mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)});
    
    avg_resp_tmp1 = mean(resp_tmp1(i, ts2use));
    avg_resp_tmp2 = mean(resp_tmp2(i, ts2use));
    
    sampling_diff(i) = (avg_resp_tmp1 - avg_resp_tmp2);%/(avg_resp_tmp1 + avg_resp_tmp2);

end
ts2use=[];
ts2use = calculate_ts2use(mean(time_stamps_v2.(groups(group2use)){1,1}), whole_window);
resp_tmp3=[];
resp_tmp3 = [resp_tmp1(:,ts2use) resp_tmp2(:,ts2use)];

sorted_idx = [];
[~,sorted_idx] = sort(sampling_diff);


rowmin = min(resp_tmp3,[],2);
rowmax = max(resp_tmp3,[],2);
Brow=[];
Brow = rescale(resp_tmp3,0,1,'InputMin',rowmin,'InputMax',rowmax);

time_stamps2use = mean(time_stamps_v2.(groups(group2use)){1,1});
colorlim2use = [0 1];

normalized_activity_taste1 = Brow(sorted_idx,1:(size(Brow,2)/2));
normalized_activity_taste2 = Brow(sorted_idx,1+(size(Brow,2)/2):end);

norm_diff = normalized_activity_taste1 - normalized_activity_taste2;
ts2use2 = calculate_ts2use(time_stamps2use(ts2use), resp_window);

% heatmaps delay responsive
figure
subplot(1,2,1)
imagesc(time_stamps2use(ts2use), 1:size(resp_tmp1,1), normalized_activity_taste1, colorlim2use)
hold on
ylim2use = [0 size(resp_tmp2,1)];
line([0 0],ylim2use, 'linestyle','--','color','k')
cb = colorbar('TickDirection','out','Ticks',[-1 0 1]);
ylabel(cb, 'Normalized \DeltaF/f');
box off
% axis square
title('Sucrose Trials')
ylabel('Neurons')
xlabel('Time (sec)')
xlim(whole_window)
yticks(round(size(resp_tmp2,1),-1))
ylim(ylim2use)

subplot(1,2,2)
imagesc(time_stamps2use(ts2use), 1:size(resp_tmp2,1),normalized_activity_taste2 , colorlim2use)
hold on
ylim2use = [0 size(resp_tmp2,1)];
line([0 0],ylim2use, 'linestyle','--','color','k')
cb = colorbar('TickDirection','out','Ticks',[-1 0 1]);
ylabel(cb, 'Normalized \DeltaF/f');
box off
% axis square
title('NaCl Trials')
ylabel('Neurons')
xlabel('Time (sec)')
xlim(whole_window)
ylim(ylim2use)
yticks(round(size(resp_tmp2,1),-1))
sgtitle('Delay Responsive Neurons','fontsize',10)

