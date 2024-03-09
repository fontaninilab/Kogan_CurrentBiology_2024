%% setup data - correct vs error
%clear all
clear; clc; close all
%add current folder and subfolders containing utility functions and data
%files
addpath(genpath(pwd))
load('neural_data_psychometrics_correct_error.mat')
groups = string(fieldnames(neural_data_v2));
groups2 = string(fieldnames(trialID_v2));

%% find responsive, correct trials only
baseline_window = [-4 -2.5; -6.5 -5; -4 -2.5; -6.5 -5]; 
resp_window_groups = [0 1.5; -1.5 0; 0 1.5; -1.5 0];

%psychometric data
[responsive_neurons, responsive_neurons_indv_tastes, pvals_all] = ...
    find_responsive_neurons_correct_error(neural_data_v2, time_stamps_v2, trialID_v2, baseline_window, resp_window_groups);

%% find selective, correct trials only
resp_window_groups = [0 1.5; -1.5 0; 0 1.5; -1.5 0];
%psychometric data
[stimulus_selective_neurons, pvals_selective] = ...
    find_selective_neurons_correct_error(neural_data_v2, time_stamps_v2,trialID_v2, resp_window_groups);

%%
resp_window_groups = repmat([0 1.5; -1 0],2,1);

%psychometric sessions
error_selective_neurons_same_stimulus=[]; error_selective_neurons_same_choice=[];
selective_neurons_same_stimulus_by_direction=[]; selective_neurons_same_choice_by_direction=[];

[error_selective_neurons_same_stimulus, error_selective_neurons_same_choice, ...
    selective_neurons_same_stimulus_by_direction, selective_neurons_same_choice_by_direction] = ...
    find_selective_neurons_correct_error_v2(neural_data_v2, time_stamps_v2, trialID_v2, resp_window_groups);

groups = string(fieldnames(neural_data_v2));
for p=1:4
group2use=p;
for h=1:2
side2use = h;
prop_error_selective(h,p) = length(intersect(responsive_neurons_indv_tastes.(groups(group2use)){side2use}, selective_neurons_same_stimulus_by_direction.(groups(group2use)){side2use}))/...
    length(responsive_neurons_indv_tastes.(groups(group2use)){side2use});
end
end
[prop_error_selective; mean(prop_error_selective)]




%% plot figure 3 - correct vs error pre learning only
% line2 = @(b,x) (b(1)*x) + b(2); %linear function
line2 = @(b,x) (b(1)*x); %linear function
slopes2 = []; r2 = [];
% selectivity scatter plots from svoboda papers
selectivity_numbers_for_stats ={}; r_sq=[];
resp_window_groups = repmat([0 1.5; -1 0],2,1);
% titles2use = ["Sampling - Pre","Delay - Pre","Sampling - Post","Delay - Post"];
titles2use = ["Sampling","Delay"];

figure
for group2use = 1:2
    
    neurons2use = []; selectivity_correct=[]; selectivity_error=[];
%         neurons2use = intersect(responsive_neurons.(groups(group2use)), stimulus_selective_neurons.(groups(group2use)){3});
%         neurons2use = stimulus_selective_neurons.(groups(group2use)){3};
    
            neurons2use = responsive_neurons.(groups(group2use));
%     neurons2use = 1:length(neural_data_v2.(groups(group2use)));
    for i=1:length(neurons2use)
        ts2use = calculate_ts2use(mean(time_stamps_v2.(groups(group2use)){1,neurons2use(i)}), resp_window_groups(group2use,:));
        
        correct_trials_suc=[]; correct_trials_NaCl=[]; error_trials_suc=[]; error_trials_NaCl=[];
        
        correct_trials_suc = find(trialID_v2.(groups(group2use)){1,neurons2use(i)}==1);
        correct_trials_NaCl = find(trialID_v2.(groups(group2use)){2,neurons2use(i)}==1);
        
        error_trials_suc = find(trialID_v2.(groups(group2use)){1,neurons2use(i)}==0);
        error_trials_NaCl = find(trialID_v2.(groups(group2use)){2,neurons2use(i)}==0);
        
        
        suc_correct_activity = mean(mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(correct_trials_suc, ts2use),2));
        NaCl_correct_activity = mean(mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(correct_trials_NaCl, ts2use),2));
        selectivity_correct(i) = (suc_correct_activity - NaCl_correct_activity);

        
        suc_error_activity = mean(mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(error_trials_suc, ts2use),2));
        NaCl_error_activity = mean(mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(error_trials_NaCl, ts2use),2));

        selectivity_error(i) =  mean([(suc_correct_activity - suc_error_activity);(NaCl_error_activity - NaCl_correct_activity)]);
        
        suc_correct_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(correct_trials_suc, :));
        NaCl_correct_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(correct_trials_NaCl, :));
        if length(error_trials_suc)>1
            suc_error_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(error_trials_suc, :));
        else
            suc_error_activity_timecourse{group2use}(i,:) = (neural_data_v2.(groups(group2use)){1,neurons2use(i)}(error_trials_suc, :));
        end
        if length(error_trials_NaCl)>1
            NaCl_error_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(error_trials_NaCl, :));
        else
            NaCl_error_activity_timecourse{group2use}(i,:) = (neural_data_v2.(groups(group2use)){2,neurons2use(i)}(error_trials_NaCl, :));
        end


    end
    
    fit_vals = nlinfit( selectivity_correct, selectivity_error,line2, [ 1 0 ]);

%     fit_vals = nlinfit( selectivity_error, selectivity_correct,line2, [ 1 0]);
    
    plot_cutoff = .5;
    x2_fits = -plot_cutoff:.01:plot_cutoff;
    y2_vals = line2(fit_vals, x2_fits);
    plts2use = [5 6];
    subplot(3,2,plts2use(group2use))
    hold on
    scatter(selectivity_correct, selectivity_error, 6,'filled','k')
    if group2use==1
        neuron2use_sampling = 432;  %14, 15, 33, 44, 46, 58  ;nonselective - 30 , alternatives - 432 (best), 305, 64, 61 (less good); 103(choice selective)
    scatter(selectivity_correct(neuron2use_sampling), selectivity_error(neuron2use_sampling), 12,'filled','c')
    %data for sampling only
    sample_selective_correct = selectivity_correct;
    sample_selective_error = selectivity_error;
    end
    if group2use==2
        %         neuron2use_delay = find(resp_delay_idx); %plot all sig fig neurons
        neuron2use_delay = 41; %plot example sig fit neuron %41 (currently used)
        scatter(selectivity_correct(neuron2use_delay), selectivity_error(neuron2use_delay), 12,'filled','c')
        delay_selective_correct = selectivity_correct;
        delay_selective_error = selectivity_error;
    end
%     7, 8, 18, 20, 21, 26, 27, 29, 30, 32, 33, 41, 42, 51

%     scatter(selectivity_error, selectivity_correct, 6,'filled','k')

    plot(x2_fits, y2_vals)
    line([-plot_cutoff plot_cutoff], [-plot_cutoff plot_cutoff], 'color','k','linestyle','--')
    axis square
    box off
    grid on
    title(titles2use(group2use))
    ylim([-plot_cutoff plot_cutoff])
    xlim([-plot_cutoff plot_cutoff])
    xlabel('Correct Selectivity')
    ylabel('Error Selectivity')

    [R,P, RL, RU] = corrcoef(  selectivity_correct, selectivity_error);
%     [R,P, RL, RU] = corrcoef( selectivity_error, selectivity_correct);

    r_sq(1, group2use) = R(1,2);
%     r_sq(2, group2use) = RL(1,2);
%     r_sq(3, group2use) = RU(1,2);
    r_sq(2, group2use) = P(1,2);
    P(1,2)
    slopes(group2use) = fit_vals(1);
    selectivity_correct_avg(group2use) = mean(abs(selectivity_correct));
    
    selectivity_numbers_for_stats{group2use} = abs(selectivity_correct);
    
%     [min(selectivity_correct) max(selectivity_correct) min(selectivity_error) max(selectivity_error)]
    mdl = fitlm(selectivity_correct, selectivity_error);
    slopes2(group2use) = mdl.Coefficients.Estimate(2)
    r2(group2use) = mdl.Rsquared.Ordinary
%     pval(group2use) = mdl.Coefficients.pValue(2)
end
% r_sq
% slopes
% [~,pval] = ttest2(selectivity_numbers_for_stats{2},selectivity_numbers_for_stats{4})


%
%plot example neurons - sampling correct vs error
neuron2use = neuron2use_sampling; %good neurons - 14, 15, 33, 44, 46, 58  ;nonselective - 30 
group2use = 1;
timestamps2use = mean(time_stamps_v2.(groups(group2use)){1,neuron2use});
plts2use = [1 3];
% figure
for i=1:2
    subplot(3,2,plts2use(i))
    hold on
    if i==1
        plot(timestamps2use, suc_correct_activity_timecourse{group2use}(neuron2use,:),'linewidth',2,'color','c')
        plot(timestamps2use, NaCl_correct_activity_timecourse{group2use}(neuron2use,:),'linewidth',2,'color','k')
        % ylim2use = get(gca, 'ylim');
        title('Correct Trials')
        
    else
        plot(timestamps2use, suc_error_activity_timecourse{group2use}(neuron2use,:),'linewidth',2,'color','c','linestyle','--')
        plot(timestamps2use, NaCl_error_activity_timecourse{group2use}(neuron2use,:),'linewidth',2,'color','k','linestyle','--')
        title('Error Trials')
        
    end
    xlim([-1 3])
    ylim2use = get(gca, 'ylim');
    line([0 0],ylim2use,'linestyle','--','color','k')
    legend('Suc','NaCl','location','northeast')

    box off
    axis square
    
    ylim(ylim2use)
    % title('Example Neuron')
    ylabel('\DeltaF/f')
    xlabel('Time (sec)')
end
% sgtitle(['Example Neuron', newline, 'Stimulus Selective'],'fontsize',12)

%
%plot example neurons - delay correct vs error
neuron2use = neuron2use_delay; %good neurons - 7, 8, 18, 20, 21, 26, 27, 29, 30, 32, 33, 41, 42, 51
timestamps2use = mean(time_stamps_v2.(groups(2)){1,neuron2use});
plts2use = [2 4];
% figure
for i=1:2
    subplot(3,2,plts2use(i))
    hold on
    if i==1
        plot(timestamps2use, suc_correct_activity_timecourse{2}(neuron2use,:),'linewidth',2,'color','c')
        plot(timestamps2use, NaCl_correct_activity_timecourse{2}(neuron2use,:),'linewidth',2,'color','k')
        % ylim2use = get(gca, 'ylim');
        title('Correct Trials')
%         legend('Sucrose - Correct','NaCl - Correct')
        
        
    else
        plot(timestamps2use, suc_error_activity_timecourse{2}(neuron2use,:),'linewidth',2,'color','c','linestyle','--')
        plot(timestamps2use, NaCl_error_activity_timecourse{2}(neuron2use,:),'linewidth',2,'color','k','linestyle','--')
                title('Error Trials')
%         legend('Sucrose - Error','NaCl - Error')
        
        
    end
    ylim2use = get(gca, 'ylim');
    line([0 0],ylim2use,'linestyle','--','color','k')
    legend('Suc','NaCl','location','northwest')

    box off
    axis square
    xlim([-3 1])
    ylim(ylim2use)
    % title('Example Neuron')
    ylabel('\DeltaF/f')
    xlabel('Time (sec)')
end
% sgtitle(['Example Neuron', newline, 'Choice Selective'],'fontsize',12)

%% ID sampling neurons
% [selectivity_correct' selectivity_error']
% [~, idx] = min(abs(sample_selective_correct + 0.42228))
% [~, idx] = min(abs(sample_selective_correct + 0.0947))

%% ID delay neurons
% [~, idx] = min(abs(delay_selective_correct - 0.16023))


%% plot figure  correct vs error pre vs post learning 
% line2 = @(b,x) (b(1)*x) + b(2); %linear function
line2 = @(b,x) (b(1)*x); %linear function
slopes2 = []; r2 = [];
% selectivity scatter plots from svoboda papers
selectivity_numbers_for_stats ={}; r_sq=[];
selectivity_correct_avg = []; selectivity_correct_sem = []; 
selectivity_error_avg = [];selectivity_error_sem = [];
    
resp_window_groups = repmat([0 1.5; -1 0],2,1);
titles2use = ["Sampling - Pre","Delay - Pre","Sampling - Post","Delay - Post"];
% titles2use = ["Sampling","Delay"];
    suc_correct_activity_timecourse=[]; NaCl_correct_activity_timecourse=[];
    suc_error_activity_timecourse=[]; NaCl_error_activity_timecourse=[];
figure
for group2use = 1:4
    
    neurons2use = []; selectivity_correct=[]; selectivity_error=[]; 

%         neurons2use = intersect(responsive_neurons.(groups(group2use)), stimulus_selective_neurons.(groups(group2use)){3});
%         neurons2use = stimulus_selective_neurons.(groups(group2use)){3};
    
            neurons2use = responsive_neurons.(groups(group2use));
%     neurons2use = 1:length(neural_data_v2.(groups(group2use)));
    for i=1:length(neurons2use)
        ts2use = calculate_ts2use(mean(time_stamps_v2.(groups(group2use)){1,neurons2use(i)}), resp_window_groups(group2use,:));
        
        correct_trials_suc=[]; correct_trials_NaCl=[]; error_trials_suc=[]; error_trials_NaCl=[];
        
        correct_trials_suc = find(trialID_v2.(groups(group2use)){1,neurons2use(i)}==1);
        correct_trials_NaCl = find(trialID_v2.(groups(group2use)){2,neurons2use(i)}==1);
        
        error_trials_suc = find(trialID_v2.(groups(group2use)){1,neurons2use(i)}==0);
        error_trials_NaCl = find(trialID_v2.(groups(group2use)){2,neurons2use(i)}==0);
        
        
        suc_correct_activity = mean(mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(correct_trials_suc, ts2use),2));
        NaCl_correct_activity = mean(mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(correct_trials_NaCl, ts2use),2));
        
        selectivity_correct(i) = (suc_correct_activity - NaCl_correct_activity);
        
        suc_error_activity = mean(mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(error_trials_suc, ts2use),2));
        NaCl_error_activity = mean(mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(error_trials_NaCl, ts2use),2));
%         
%         if selectivity_correct(i)>=0
%             selectivity_error(i) =  (suc_correct_activity - suc_error_activity);
%         else
%             selectivity_error(i) =  (NaCl_error_activity - NaCl_correct_activity);
%             
%         end
        selectivity_error(i) =  mean([(suc_correct_activity - suc_error_activity);(NaCl_error_activity - NaCl_correct_activity)]);
        
        suc_correct_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(correct_trials_suc, :));
        NaCl_correct_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(correct_trials_NaCl, :));
        if length(error_trials_suc)>1
            suc_error_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){1,neurons2use(i)}(error_trials_suc, :));
        else
            suc_error_activity_timecourse{group2use}(i,:) = (neural_data_v2.(groups(group2use)){1,neurons2use(i)}(error_trials_suc, :));
        end
        if length(error_trials_NaCl)>1
            NaCl_error_activity_timecourse{group2use}(i,:) = mean(neural_data_v2.(groups(group2use)){2,neurons2use(i)}(error_trials_NaCl, :));
        else
            NaCl_error_activity_timecourse{group2use}(i,:) = (neural_data_v2.(groups(group2use)){2,neurons2use(i)}(error_trials_NaCl, :));
        end


    end
    selectivity_all{group2use} = [selectivity_correct; selectivity_error]';

    
    fit_vals = nlinfit( selectivity_correct, selectivity_error,line2, [ 1 0 ]);

    
    
    
    plot_cutoff = .5;
    x2_fits = -plot_cutoff:.01:plot_cutoff;
    y2_vals = line2(fit_vals, x2_fits);
    plts2use = [5 6];
    subplot(2,2,group2use)
    hold on
    scatter(selectivity_correct, selectivity_error, 6,'filled','k')
%     scatter(selectivity_error, selectivity_correct, 6,'filled','k')

    plot(x2_fits, y2_vals)
    line([-plot_cutoff plot_cutoff], [-plot_cutoff plot_cutoff], 'color','k','linestyle','--')
    axis square
    box off
    grid on
    title(titles2use(group2use))
    ylim([-plot_cutoff plot_cutoff])
    xlim([-plot_cutoff plot_cutoff])
    xlabel('Correct Selectivity')
    ylabel('Error Selectivity')

    [R,P, RL, RU] = corrcoef(  selectivity_correct, selectivity_error);

    r_sq(1, group2use) = R(1,2);

    r_sq(2, group2use) = P(1,2);
    P(1,2)
    slopes(group2use) = fit_vals(1);
    selectivity_correct_avg(group2use) = mean(abs(selectivity_correct));
    selectivity_correct_sem(group2use) = find_sem(abs(selectivity_correct));

    selectivity_error_avg(group2use) = mean(abs(selectivity_error));
    selectivity_error_sem(group2use) = find_sem(abs(selectivity_error));

    selectivity_numbers_for_stats{group2use} = abs(selectivity_correct);
    selectivity_numbers_for_stats_error{group2use} = abs(selectivity_error);

    mdl = fitlm(selectivity_correct, selectivity_error);
    slopes2(group2use) = mdl.Coefficients.Estimate(2);
    r2(group2use) = mdl.Rsquared.Ordinary;
end
slopes2
% selectivity_correct_avg
%% testing significant difference in correlation between sampling and delay
num_resamples = 10000;

%pairs of groups - 
% row 1 = prelearning sample vs delay
% row 2 = postlearning sample vs delay
% row 3 = sample prelearning vs postlearning
% row 4 = delay prelearning vs postlearning
titles_for_comparisons = ["Prelearning\nSample vs Delay", "Postlearning\nSample vs Delay",...
    "Sample\nPrelearning vs Postlearning","Delay\nPrelearning vs Postlearning"];
pairs2test = [1 2; 3 4;1 3; 2 4 ]; 

corr_diff=zeros(size(pairs2test,1),num_resamples);
for p=1:size(pairs2test,1)
    num_datapoints2use = min([size(selectivity_all{pairs2test(p,1)},1) size(selectivity_all{pairs2test(p,2)},1)]);
    for i=1:num_resamples
        data_null = [selectivity_all{pairs2test(p,1)}; selectivity_all{pairs2test(p,2)}];
        X1 = datasample(data_null, num_datapoints2use);
        r1 = corr(X1(:,1),X1(:,2))^2;
        X2 = datasample(data_null , num_datapoints2use);
        r2 = corr(X2(:,1),X2(:,2))^2;

        corr_diff(p,i) = (r2- r1);
    end
end
%%
figure
tiledlayout("flow","TileSpacing","compact")
actual_corr_diff=zeros(1,size(pairs2test,1));
for p=1:size(pairs2test,1)
    nexttile
    histogram(corr_diff(p,:), 'Normalization','probability')
    hold on
    ylim2use = get(gca, 'ylim');
    X1_all = selectivity_all{pairs2test(p,1)};
    r1_all = corr(X1_all(:,1),X1_all(:,2))^2;
    X2_all = selectivity_all{pairs2test(p,2)};
    r2_all = corr(X2_all(:,1),X2_all(:,2))^2;

    actual_corr_diff(p) =  (r2_all - r1_all);
    line([actual_corr_diff(p) actual_corr_diff(p)], ylim2use)
    title(compose(titles_for_comparisons(p)))
    xlim([-.3 .5])
    ylabel('Percentage')

    pval_distribution(p) = (length(find(corr_diff(p,:)>actual_corr_diff(p)))/num_resamples);
end
% nexttile
% legend("Actual")
% pd = fitdist(corr_difference', 'Normal');
% ci = paramci(pd)
% length(find(corr_difference>actual_r2))
%%
%correct stats
[~,pval] = ttest2(selectivity_numbers_for_stats{1},selectivity_numbers_for_stats{3})
[~,pval] = ttest2(selectivity_numbers_for_stats{2},selectivity_numbers_for_stats{4})

%error stats
[~,pval] = ttest2(selectivity_numbers_for_stats_error{1},selectivity_numbers_for_stats_error{3})
[~,pval] = ttest2(selectivity_numbers_for_stats_error{2},selectivity_numbers_for_stats_error{4})

%% mean selectivity numbers
for i=1:4
selectivity_numbers_avg(i) = mean(selectivity_numbers_for_stats{i})
end
%% plot mean selectivity - pre vs post learning - delay
figure
subplot(1,2,1)
xvals2 = bar_custom_mixtures(reshape(selectivity_correct_avg,2,[]));
hold on
errorbar(xvals2', reshape(selectivity_correct_avg,2,[]),[0 0; 0 0], reshape(selectivity_correct_sem,2,[]),'linestyle','none','color','k')
box off
axis square
xticklabels(["Sampling","Delay"])
title('Correct Selectivity')
ylabel('\Delta\DeltaF/f')
ylim([0 .08])
yticks([0 .04 .08])


grps2use = [2 4];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

subplot(1,2,2)
hold on
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})),'linewidth',2,'color','k')
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})),'linewidth',2,'color','c')
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})))
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})))

xlim([-3 1])
ylim2use = [0.02 .11];
% ylim2use = get(gca, 'ylim');
line([0 0], ylim2use, 'linestyle','--','color','k')
legend('Pre','Post','location','northwest')
box off
axis square
ylabel('\Delta\DeltaF/f')
xlabel('Time (sec)')
title('Delay')
ylim(ylim2use)
%% sig difference between lines, timecourse - delay
grps2use = [2 4];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

pre_learning_data = abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)});
post_learning_data = abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)});
p=[];
ts2use = calculate_ts2use(timestamps2use, [-3 1]);
for i=1:length(ts2use)
    [~,p(i)] = ttest2(pre_learning_data(:,ts2use(i)), post_learning_data(:,ts2use(i)));
end

timestamps2use(ts2use(find(p<0.05/length(ts2use))))

%% plot mean selectivity - pre vs post learning - sampling
figure
subplot(1,2,1)
xvals2 = bar_custom_mixtures(reshape(selectivity_correct_avg,2,[]));
hold on
errorbar(xvals2', reshape(selectivity_correct_avg,2,[]),[0 0; 0 0], reshape(selectivity_correct_sem,2,[]),'linestyle','none','color','k')
box off
axis square
xticklabels(["Sampling","Delay"])
title('Correct Selectivity')
ylabel('\Delta\DeltaF/f')
ylim([0 .08])
yticks([0 .04 .08])

subplot(1,2,2)
grps2use = [1 3];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

subplot(1,2,2)
hold on
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})),'linewidth',2,'color','k')
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})),'linewidth',2,'color','c')
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})))
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})))

xlim([-1 3])
ylim2use = [0.02 .11];
% ylim2use = get(gca, 'ylim');
line([0 0], ylim2use, 'linestyle','--','color','k')
legend('Pre','Post','location','northwest')
box off
axis square
ylabel('\Delta\DeltaF/f')
xlabel('Time (sec)')
title('Sampling')
ylim(ylim2use)
%% sig difference between lines, timecourse - sampling
grps2use = [1 3];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

pre_learning_data = abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)});
post_learning_data = abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)});
p=[];
ts2use = calculate_ts2use(timestamps2use, [-1 3]);
for i=1:length(ts2use)
    [~,p(i)] = ttest2(pre_learning_data(:,ts2use(i)), post_learning_data(:,ts2use(i)));
end

timestamps2use(ts2use(find(p<0.05/length(ts2use))))
%% time to half maximum for delay selectivity
grps2use = [2 4];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

pre_selectivity_timecourse = mean(abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)}));
post_selectivity_timecourse = mean(abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)}));

%max selectivity in delay
delay_window = calculate_ts2use(timestamps2use, [-1.5, 0]);
[~,max_idx_pre] = max(pre_selectivity_timecourse(delay_window));
[~,max_idx_post] = max(post_selectivity_timecourse(delay_window));

%mean selectivity in delay, or mid point of selectivity across delay
half_max_pre = (pre_selectivity_timecourse(delay_window(max_idx_pre))+min(pre_selectivity_timecourse(delay_window)))/2;
half_max_post = (post_selectivity_timecourse(delay_window(max_idx_post))+min(post_selectivity_timecourse(delay_window)))/2;

%index of time stamp where selectivity exceeds half max
half_max_idx_pre = find(pre_selectivity_timecourse(delay_window)>=half_max_pre);
half_max_idx_post = find(post_selectivity_timecourse(delay_window)>=half_max_post);

%time point in delay where selectivity exeeds half max selectivity
half_max_timestamp_pre = timestamps2use(delay_window(half_max_idx_pre(1)))
half_max_timestamp_post = timestamps2use(delay_window(half_max_idx_post(1)))

%% time to half maximum for delay selectivity - stats with bootstrap approach
grps2use = [2 4];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});
delay_window = calculate_ts2use(timestamps2use, [-1.5, 0]);

pre_selectivity_timecourse = (abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)}));
post_selectivity_timecourse = (abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)}));

num_observations2use = min([size(pre_selectivity_timecourse,1) size(post_selectivity_timecourse,1)]);
data_null_selectivity_timecourse = [pre_selectivity_timecourse; post_selectivity_timecourse];

num_resamples = 10000; 
difference_in_halfmax_bootstrapped = zeros(1,num_resamples);
for i=1:num_resamples
    %subsampled populations timecourses
    subsample1 = randperm(size(data_null_selectivity_timecourse,1), size(pre_selectivity_timecourse,1));
    subsample2 = setdiff(1:size(data_null_selectivity_timecourse,1), subsample1);

    X1_resample = mean(data_null_selectivity_timecourse(subsample1,:));
    X2_resample = mean(data_null_selectivity_timecourse(subsample2,:));

    % % X1_resample = mean(datasample(data_null_selectivity_timecourse, num_observations2use,1,'Replace',false));
    % X2_resample = mean(datasample(data_null_selectivity_timecourse, num_observations2use,1));


    %max selectivity in delay
    [~,max_idx_X1] = max(X1_resample(delay_window));
    [~,max_idx_X2] = max(X2_resample(delay_window));

    %mean selectivity in delay, or mid point of selectivity across delay
    half_max_X1 = (X1_resample(delay_window(max_idx_X1))+min(X1_resample(delay_window)))/2;
    half_max_X2 = (X2_resample(delay_window(max_idx_X2))+min(X2_resample(delay_window)))/2;

    %index of time stamp where selectivity exceeds half max
    half_max_idx_X1 = find(X1_resample(delay_window)>=half_max_X1);
    half_max_idx_X2 = find(X2_resample(delay_window)>=half_max_X2);

    %time point in delay where selectivity exeeds half max selectivity
    half_max_timestamp_X1 = timestamps2use(delay_window(half_max_idx_X1(1)));
    half_max_timestamp_X2 = timestamps2use(delay_window(half_max_idx_X2(1)));

    %difference in half max between subsampled populations
    difference_in_halfmax_bootstrapped(i) = abs(half_max_timestamp_X2 - half_max_timestamp_X1);
end
actual_difference_halfmax_selectivity = abs(half_max_timestamp_pre - half_max_timestamp_post);
figure
histogram(difference_in_halfmax_bootstrapped,'Normalization','probability')
ylim2use = get(gca,'ylim');
line([actual_difference_halfmax_selectivity actual_difference_halfmax_selectivity], ylim2use, 'linestyle',':','color','r')

% figure
% bar([mean(difference_in_halfmax_bootstrapped) actual_difference_halfmax_selectivity])
%%
figure
hold on
plot(timestamps2use, X1_resample)
plot(timestamps2use, X2_resample)



%%
grps2use = [1 3];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

figure
subplot(1,3,2)

hold on
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})),'linewidth',2,'color','k')
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})),'linewidth',2,'color','c')
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})))
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})))

xlim([-1 2])
ylim2use = [0.02 .05];
% ylim2use = get(gca, 'ylim');
line([0 0], ylim2use, 'linestyle','--','color','k')
legend('Pre','Post','location','northwest')
box off
axis square
ylabel('\Delta\DeltaF/f')
xlabel('Time (sec)')
title('Sampling')
ylim(ylim2use)
%
grps2use = [2 4];
timestamps2use = mean(time_stamps_v2.(groups(grps2use(1))){1,1});

subplot(1,3,3)
hold on
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})),'linewidth',2,'color','k')
plot(timestamps2use, mean(abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})),'linewidth',2,'color','c')
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(1)} - NaCl_correct_activity_timecourse{grps2use(1)})))
plot_sem(timestamps2use, (abs(suc_correct_activity_timecourse{grps2use(2)} - NaCl_correct_activity_timecourse{grps2use(2)})))

xlim([-2 1])
ylim2use = [0.02 .1];
% ylim2use = get(gca, 'ylim');
line([0 0], ylim2use, 'linestyle','--','color','k')
legend('Pre','Post','location','northwest')
box off
axis square
ylabel('\Delta\DeltaF/f')
xlabel('Time (sec)')
title('Delay')
ylim(ylim2use)

subplot(1,3,1)
xvals2 = bar_custom_mixtures(reshape(selectivity_correct_avg,2,[]));
hold on
errorbar(xvals2', reshape(selectivity_correct_avg,2,[]),[0 0; 0 0], reshape(selectivity_correct_sem,2,[]),'linestyle','none','color','k')
box off
axis square
xticklabels(["Sampling","Delay"])
title('Selectivity - Responsive Neurons')
ylabel('\Delta\DeltaF/f')
ylim([0 .08])
yticks([0 .04 .08])