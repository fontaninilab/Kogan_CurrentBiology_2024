%% curve fitting on neural data - Figs 3 and 7
% cd('D:\Code_data_analysis\Curr_bio_data_upload')
%clear all
clear; clc; close all
%add current folder and subfolders containing utility functions and data
%files
addpath(genpath(pwd))
load('Figure_3_7_neural_data_mixtures_v3.mat')
%
groups = string(fieldnames(neural_data));
%% define functions
line2 = {@(b,x) (b(1)*x./100) + b(2),... %linear function
            @(b,x) b(1)+(b(2)-b(1))./(1+10.^((50-x)*b(3))),...%sig function
            @(b,x) (b(1)*abs(b(1)*x-b(1)*50)/100)+b(2)}; %V shape
%% run fits on all responsive neurons
%this takes a long time to run, output is also saved and used for plotting
X2 = 100:-1:0;
resp_window_groups = [0 1.5; -1.5 0; 0 1.5; -1.5 0];
% bins2use = -1:.5:2;
i=1;
bins2use_groups = [-1:.5:3.5; -3:.5:1.5; -1:.5:3.5; -3:.5:1.5];

for group2use=1:4
    for p=1:size(bins2use_groups,2)-1
        resp_window2use = [bins2use_groups(group2use,p) bins2use_groups(group2use,p+1)];
        data2fit.(groups(group2use)){p} = avg_data_curve_fitting_mixtures_v2(neural_data, time_stamps, resp_window2use, group2use);
        
        [SS_resid.(groups(group2use)), ~, fit_vals.(groups(group2use)){p}] = calculate_fits_mixtures(data2fit.(groups(group2use)){p});
        
        [sig_neurons_f.(groups(group2use)){p}, f_stat.(groups(group2use)){p}] = calculate_f_stat_mixtures_fit(SS_resid.(groups(group2use)));
    end
end
%% identify responsive neurons by six tastes - use for example neurons only
baseline_window = [-4 -2.5; -6.5 -5; -4 -2.5; -6.5 -5]; 
resp_window_groups = [0 1.5; -1.5 0; 0 1.5; -1.5 0];
[responsive_neurons, responsive_neurons_indv_tastes, pvals_all] = find_responsive_neurons_v2(neural_data, time_stamps, baseline_window, resp_window_groups);
%
%% save curve fitting results
%they are loaded in below
% save('Neural_data_curve_fitting_results.mat','data2fit','SS_resid','fit_vals','sig_neurons_f','f_stat','responsive_neurons')
%% load saved curve fitting results
clear; clc; close all

% cd('D:\Code_data_analysis\Curr_bio_data_upload')
addpath(genpath(pwd))
load('Neural_data_curve_fitting_results.mat')
load('Figure_2_neural_data_w_significant_neurons.mat')
load('proportions_resp_selective.mat')
load('Figure_3_7_neural_data_mixtures_v3.mat')
groups = string(fieldnames(neural_data));
bins2use_groups = [-1:.5:3.5; -3:.5:1.5; -1:.5:3.5; -3:.5:1.5];
resp_window_groups = [0 1.5; -1.5 0; 0 1.5; -1.5 0];

%% curve fitting figure - prelearning
% %plot examples
neurons2use=[];
group2use=1;
bin2use=6;
% bins2use_groups(group2use, bin2use)
neurons2use = intersect(responsive_neurons.(groups(group2use)), sig_neurons_f.(groups(group2use)){bin2use}{2,1});
% neuron2use = neurons2use(27); %10, 11, 14, 16, 27*
neuron2use = 2227; %2227;%953;%;%1139;
window2use_plt = [-1 3];
subplts2use=[1 2];
figure
plot_example_fits(neural_data, data2fit, f_stat, fit_vals, time_stamps, neuron2use,  group2use, window2use_plt, [0 1.5], bin2use, subplts2use)
% plot_example_fits(neural_data, data2fit, f_stat, fit_vals, time_stamps,  neuron2use,  group2use, window2use_plt, resp_window)
% %plot example sigmoid
group2use=2;
bin2use=6;
% bins2use_groups(group2use, bin2use)
neurons2use = intersect(responsive_neurons.(groups(group2use)), sig_neurons_f.(groups(group2use)){bin2use}{2,2});
neuron2use = neurons2use(14);%14 %16
% neuron2use = 13;
window2use_plt = [-3 1];
subplts2use=[3 4];
% figure
plot_example_fits(neural_data, data2fit, f_stat, fit_vals, time_stamps, neuron2use,  group2use, window2use_plt, [-1.5 0], bin2use, subplts2use)

% proportion timecourse
groups2use = [1 2];
for h=1:length(groups2use)
    for bin2use = 1:size(bins2use_groups,2)-1
        for j=1:3
            sig_fits_bin = length(intersect(responsive_neurons_v2.(groups(h)), sig_neurons_f.(groups(h)){bin2use}{2,j}));
            responsive_neurons_bin = length(responsive_neurons_v2.(groups(h)));
            proportion_fits(j, bin2use, h) = sig_fits_bin/responsive_neurons_bin;
        end
    end
end
% plotting
titles = ["Sampling","Delay"];
order = [1 3 2 4];
% figure('Position',psn)
% figure
colors2use=[];
colors2use = brewermap(3,'paired');
time_course_order = [1 3];
for h=1:length(groups2use)
   subplot(2,4,time_course_order(h)+4)
   hold on
   x = bins2use_groups(h,1:end-1);
   y = proportion_fits(:,:,h);
   for q=1:2 %changed from 3
       plot(x,y(q,:),'color', colors2use(q,:),'linewidth',2);
   end
   box off
   if ismember(h,[1 3])
       xlim([-1 3])
       xticks(-1:3)
   else
      xlim([-3 1]) 
      xticks(-3:1)
   end
   line([0 0],[0 .25],'color','k','linestyle','--')

   if ismember(h,1)
       legend('Line','Sigmoid','Taste','Location','northeastoutside')
   
   else
       legend('Line','Sigmoid','Choice','Location','northeastoutside')
   end
   ylim([0 .22])
   yticks(0:.1:.2)
   title(titles(h))
   ylabel('Proportion')
   xlabel('Time (sec)')
   axis square

end
% plot proportion of all fits in window and pie charts

max_proportion_responsive = zeros(3,2);
maxbin_idx = zeros(3,2);
% figure
pie_chart_order = [2 4];
for h=1:length(groups2use)
   x = bins2use_groups(h,1:end-1);
   y = proportion_fits(:,:,h);
   resp_window_proportions = calculate_ts2use(x,resp_window_groups(h,:));
   [max_proportion_responsive(1:2,h), maxbin_idx(1:2,h)] = max(proportion_fits(1:2,resp_window_proportions(1:end-1),h),[],2);
   max_proportion_responsive(3,h) = 1-sum(max_proportion_responsive(1:3,h));
   
   %bin w max response
   tmp1 = []; max_bins_tmp = [];
   tmp1 = repmat(resp_window_proportions(1:end-1),2,1);
   max_bins_tmp = maxbin_idx(1:2,h);
   max_fits_bin_ts(:,h) = [x(tmp1(1, max_bins_tmp(1)));... %linear bin peak
   x(tmp1(2, max_bins_tmp(2)))] %sigmoid bin peak;

   subplot(2,4,pie_chart_order(h)+4)
   colormap(colors2use)
   pie(max_proportion_responsive(:,h))
   title(titles(h))
%    if h==1
%        legend('Line','Sig','Mixture','None')
        legend('Line','Sig','None','location','northeastoutside')
%    end

end
%% proportions timecourse
groups2use = [1 2 3 4]; proportion_fits_pre_post=[];
sig_fits_bin_all=[]; resp_neurons_bin_all=[];
for h=1:length(groups2use)
    for bin2use = 1:size(bins2use_groups,2)-1
        for j=1:3
            sig_neurons_all{j,h,bin2use} = intersect(responsive_neurons_v2.(groups(h)), sig_neurons_f.(groups(h)){bin2use}{2,j});
            sig_fits_bin = length(intersect(responsive_neurons_v2.(groups(h)), sig_neurons_f.(groups(h)){bin2use}{2,j}));
            responsive_neurons_bin = length(responsive_neurons_v2.(groups(h)));
            proportion_fits_pre_post(j, bin2use, h) = sig_fits_bin/responsive_neurons_bin;
            sig_fits_bin_all(j,bin2use,h) = sig_fits_bin;
            resp_neurons_bin_all(j,bin2use,h) = responsive_neurons_bin;
        end
    end
end
%% pull out step fit neurons in delay, prelearning for correct vs error fig
bins_idx = calculate_ts2use(bins2use_groups(2,:),[-1.5 1]);
all_sig_neurons=[];
for i=1:length(bins_idx)
    all_sig_neurons = [all_sig_neurons sig_neurons_all{2,:,i}];
    
end
sigmoid_fit_neurons_delay = unique(all_sig_neurons);
p=1;
for i=1:length(responsive_neurons_v2.pre_lateral)
     if ismember(responsive_neurons_v2.pre_lateral(i), sigmoid_fit_neurons_delay)
        resp_delay_idx(i) =1;
     else
         resp_delay_idx(i) =0;
     end
end


%% stats on curve fitting timecourse
pre_learning_numbers = [sig_fits_bin_all(2,:,2); resp_neurons_bin_all(2,:,2)]';
post_learning_numbers = [sig_fits_bin_all(2,:,4); resp_neurons_bin_all(2,:,4)]';

for i=1:size(pre_learning_numbers,1)
    pvals_prevspost_sigmoid(i) = chi2test([pre_learning_numbers(i,:); post_learning_numbers(i,:)]);
end
sig_bins = find(pvals_prevspost_sigmoid<0.05);

%numbers for paper
time_bins = bins2use_groups(2,sig_bins);
pre_values = proportion_fits_pre_post(2,sig_bins,2)*100;
post_values = proportion_fits_pre_post(2,sig_bins,4)*100;
numbers_for_paper = [time_bins; pre_learning_numbers(sig_bins,:)'; pre_values; post_learning_numbers(sig_bins,:)';post_values; pvals_prevspost_sigmoid(sig_bins)]

%% prelearning proportion fits stats
%compare sampling to delay - max proportion - chi2 test
for i=1:2 %loop through linear and sigmoid fit types
    proportion_data = [];
    proportion_data = [max_proportion_responsive(i,1)*length(responsive_neurons_v2.pre_central) length(responsive_neurons_v2.pre_central);...
        max_proportion_responsive(i,2)*length(responsive_neurons_v2.pre_lateral) length(responsive_neurons_v2.pre_lateral)]
    pval(i) = chi2test(proportion_data)
end

%% compare fstats
f_stat_all_resp=[];
for group2use = 1:4
avg_fstat_tmp=[];
% group2use=4;
for p=1:size(f_stat.(groups(group2use)),2)
    neurons2use =[];
    neurons2use = intersect(responsive_neurons_v2.(groups(group2use)), sig_neurons_f.(groups(group2use)){1,p}{2,1});
    avg_fstat_tmp{1,p} = (f_stat.(groups(group2use)){1,p}(1,neurons2use));
    avg_fstat_tmp{2,p} = (f_stat.(groups(group2use)){1,p}(2,neurons2use));

end
f_stat_all_resp.(groups(group2use)) = avg_fstat_tmp;
end

% %%
% diff_f_stat_tmp=[];
% j=6;
%     diff_f_stat_tmp(1) = mean(f_stat_all_resp.(groups(2)){1,j} - f_stat_all_resp.(groups(2)){2,j})
% 
%     diff_f_stat_tmp(2) = mean(f_stat_all_resp.(groups(4)){1,j} - f_stat_all_resp.(groups(4)){2,j})
% 
% figure
% bar(diff_f_stat_tmp)
%% pre vs post comparison and plotting (figure 7)

for h=1:length(groups)
    for bin2use = 1:size(bins2use_groups,2)-1
        for j=1:3
            sig_fits_bin = length(intersect(responsive_neurons_v2.(groups(h)), sig_neurons_f.(groups(h)){bin2use}{2,j}));
            responsive_neurons_bin = length(responsive_neurons_v2.(groups(h)));
            proportion_fits(j, bin2use, h) = sig_fits_bin/responsive_neurons_bin;
        end
    end
end
max_proportion_responsive=[];
for h=1:length(groups)
   x = bins2use_groups(h,1:end-1);
   y = proportion_fits(:,:,h);
   resp_window_proportions = calculate_ts2use(x,resp_window_groups(h,:));
   max_proportion_responsive(1:2,h) = max(proportion_fits(1:2,resp_window_proportions(1:end-1),h),[],2);
   max_proportion_responsive(3,h) = 1-sum(max_proportion_responsive(1:2,h));
end
figure
subplot(1,3,1)
bar_custom_mixtures(reshape(proportion_responsive_v2,2,[]))
xlim([.5 2.5])
xticklabels(["Sampling","Delay"])
yticks([0 .15 .3])
ylim([0 .3])
axis square
legend('Pre','Post','Location','northeastoutside')
ylabel('Fraction of All Neurons')
title('Responsive Neurons')

subplot(1,3,2)
bar_custom_mixtures(reshape(proportion_selective,2,[]))
xlim([.5 2.5])
xticklabels(["Sampling","Delay"])
yticks([0 .3 .6])
ylim([0 .6])
axis square
legend('Pre','Post','Location','northeastoutside')
ylabel('Fraction of Responsive Neurons')
title('Selective Neurons')

subplot(1,3,3)

hold on
b1 = bar([1  4 ],max_proportion_responsive(1:2,[1 2])' ,'stacked','BarWidth',.3);
b2 = bar([2 5],max_proportion_responsive(1:2,[3 4])' ,'stacked','BarWidth',.3);

b1(1).FaceColor = 'flat';
b1(2).FaceColor = 'flat';
b2(1).FaceColor = 'flat';
b2(2).FaceColor = 'flat';

colors2use=[];
colors2use = brewermap(5,'paired'); 

b1(1).CData = repmat(colors2use(2,:),2,1);
b1(2).CData = repmat(colors2use(4,:),2,1);

b2(1).CData = repmat(colors2use(2,:),2,1);
b2(2).CData = repmat(colors2use(4,:),2,1);

b2(1).EdgeColor = 'c'; 
b2(2).EdgeColor = 'c'; 

b1(1).LineWidth = 2;
b1(2).LineWidth = 2;

b2(1).LineWidth = 2;
b2(2).LineWidth = 2;

xticks([1.5 4.5])
xticklabels(["Sampling","Delay"])
yticks(0:.2:.4)
ylim([0 .42])
legend('Linear - Pre','Sigmoid - Pre','Linear - Post','Sigmoid - Post','location','northeastoutside')
axis square
box off
title('Significant Fits')
ylabel('Fraction of Responsive Neurons')
%% stats for pre vs post learning comparison of fits
num_fits=[];
for i=1:length(groups)
num_fits(:,i) = max_proportion_responsive(:,i)*length(responsive_neurons_v2.(groups(i)));
end

%linear fits - delay - pre vs post
pval_pre = chi2test([max_proportion_responsive(1,2)*length(responsive_neurons_v2.(groups(2))) length(responsive_neurons_v2.(groups(2)));...
max_proportion_responsive(1,4)*length(responsive_neurons_v2.(groups(4))) length(responsive_neurons_v2.(groups(4)))]);

%sigmoid fits - delay - pre vs post
pval_post = chi2test([max_proportion_responsive(2,2)*length(responsive_neurons_v2.(groups(2))) length(responsive_neurons_v2.(groups(2)));...
max_proportion_responsive(2,4)*length(responsive_neurons_v2.(groups(4))) length(responsive_neurons_v2.(groups(4)))]);

%% category tuning
%
%CTI pre vs Post
neural_data_avg_activity=[]; CTI_all=[];

for group2use=1:4
    window2use = resp_window_groups(group2use,:);
    neural_data_avg_activity.(groups(group2use)) = avg_data_curve_fitting_mixtures_v2(neural_data, time_stamps, window2use, group2use);
    CTI_all.(groups(group2use)) = calculate_CTI_mixtures(neural_data_avg_activity.(groups(group2use)));
end
%% plot category tuning
for group2use = 1:4
    resp_neurons_CTI(group2use) = mean(CTI_all.(groups(group2use))(responsive_neurons_v2.(groups(group2use))));
    resp_neurons_CTI_sem(group2use) = find_sem(CTI_all.(groups(group2use))(responsive_neurons_v2.(groups(group2use))));
    resp_neurons_CTI_values{group2use} = CTI_all.(groups(group2use))(responsive_neurons_v2.(groups(group2use)));
end
avg_CTI = reshape(resp_neurons_CTI,2,[]);
sem_CTI = reshape(resp_neurons_CTI_sem,2,[]);

figure
x2 = bar_custom_mixtures(avg_CTI);
hold on

errorbar(x2', avg_CTI, [0 0; 0 0],sem_CTI,'color','k','linestyle','none')
xlim([.5 2.5])
xticklabels(["Sampling","Delay"])
axis square
legend('Pre','Post','Location','northeastoutside')
ylabel('CTI')
title('Category Tuning Index')
%%
%CTI stats ttest
[~, pval] = ttest2(CTI_all.(groups(1))(responsive_neurons_v2.(groups(1))),CTI_all.(groups(3))(responsive_neurons_v2.(groups(3))))
[~, pval] = ttest2(CTI_all.(groups(2))(responsive_neurons_v2.(groups(2))),CTI_all.(groups(4))(responsive_neurons_v2.(groups(4))))