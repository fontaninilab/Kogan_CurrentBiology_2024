%% load in behavior data
%clear all
clear; clc; close all
%add current folder and subfolders containing utility functions and data
%files
addpath(genpath(pwd))
load('Figure_1_6_psychometric_beh_trials_data.mat')
load('Figure_1_6_licking_analysis.mat')
load('Figure_1_6_avg_behavior_data.mat')
Animals = vertcat(total_perf_all.animal);


%% plot psychometric data
colors=[]; ic50_fit=[];
figure
%sucrose concentrations
x=[100 85 75 65 60 40 35 25 15 0];

for p=1:2
    psychometric_avg=[];
    if p==2
        psychometric=[];
        psychometric=after;
    else
        psychometric = [];
        psychometric=before;
    end
    colors = [0 0 0; 0 1 1];

    for i=1:size(psychometric,1)
        minimum = psychometric(end);
        maximum=psychometric(1);
        %fit using a sigmoid function
        sigfunc = @(A,x) minimum+(maximum-minimum)./(1+10.^((A(1)-x)*A(2))); %fix min and max values
        A_init=[ 50 0 ];  % "x50", "slope"
        sigfunc(A_init,x);
        A_fit = nlinfit(x, psychometric(i,:), sigfunc, A_init);

        x_eval = (100:-5:0);
        A_eval = sigfunc(A_fit,x_eval);
        if i==size(psychometric,1) %plot mean
            scatter(x,psychometric(i,:), 50,colors(p,:),'filled')
            ylim([0 1])
            hold on
            plot(x_eval,A_eval,'Color',colors(p,:),'LineWidth',3)
            ic50_fit(p,:) = A_fit;
            %
        end
    end

    title('Psychometric Performance')
    ylabel('Sucrose Choice')
    xlabel('Percent Sucrose')
    legend('PreTraining', 'Pretratining-Fit','PostTraining','PostTraining - Fit','Location','Northwest')
end
ic50_fit % fit values pre and post

%% plot mixture pair performance and anova, pre-learning

before3 = before;
after3 = after;

before3(:,6:end,:) = 1-before3(:,6:end);
after3(:,6:end,:) = 1-after3(:,6:end);
y11=[]; y22=[];
x = [60 65 75 85 100];
% y1 = flip(mean([before3(:,1:5);flip(before3(:,6:end),2)]));
y11(:,:,1) = before3(:,1:5);
y11(:,:,2) = flip(before3(:,6:end),2);

y1=flip(mean((mean(y11,3))));
y33 = flip(((mean(y11,3))))
% y2 = flip(mean([after3(:,1:5);flip(after3(:,6:end),2)]));
y22(:,:,1) = after3(:,1:5);
y22(:,:,2) = flip(after3(:,6:end),2);

y2=flip(mean((mean(y22,3))));
y3=flip(((mean(y22,3))))
% figure(222)
figure
scatter(x,y1,'filled')
hold on
scatter(x,y2,'filled')
ylim([.5 1])
pre_mixture_pairs = mean(y11,3);
anova1(pre_mixture_pairs)
%% calculate sampling duration and directional latency for each taste
overall_sampling_duration=[]; overall_lateral_latency=[];
avg_sampling_duration_ANOVA = []; avg_lateral_latency_ANOVA = []; taste_lbls_ANOVA = [];
avg_LR_latency = []; avg_LR_latency_ANOVA=[]; LR_lbls_ANOVA = []; LR_for_ttest=[];

for h=1:2
    avg_sampling_duration=[]; avg_lateral_latency=[];
    for z=1:6
        tastes_all_trials="";
        trial2use = before_trials(h,z).trials; %use data from prelearning
        % trial2use = before_trials(h,z).trials; %post-learning

        for i=1:length(trial2use)
            tastes_all_trials(i) = string(trial2use(i).TasteID);
        end

        tastes2use = fliplr(unique(tastes_all_trials));
        for j=1:length(tastes2use)
            p=1; lateral_licks_latency=[]; sampling_duration=[];
            for i=1:length(trial2use)
                if string(trial2use(i).TasteID) == tastes2use(j)
                    sampling_duration(p) = trial2use(i).centSp(end);
                    lateral_licks_latency(p) = min([trial2use(i).LeftSp trial2use(i).RightSp]);
                    p=p+1;
                end
            end
            avg_sampling_duration(z,j) = mean(sampling_duration);
            avg_lateral_latency(z,j) = mean(lateral_licks_latency);
        end
        tt=1; rr=1;
        for i=1:length(trial2use)

            if trial2use(i).L_R_trial == 1
                lateral_licks_latency_L(tt) = min([trial2use(i).LeftSp trial2use(i).RightSp]);
                tt=tt+1;
            end
            if trial2use(i).L_R_trial == 2
                lateral_licks_latency_R(rr) = min([trial2use(i).LeftSp trial2use(i).RightSp]);
                rr=rr+1;
            end
        end
        avg_LR_latency(z,:) = [mean(lateral_licks_latency_L) mean(lateral_licks_latency_R) ]
    end
    overall_sampling_duration(h,:) = mean(avg_sampling_duration);
    overall_lateral_latency(h,:) = mean(avg_lateral_latency);

    avg_sampling_duration_ANOVA = [avg_sampling_duration_ANOVA; reshape(avg_sampling_duration,[],1)];
    avg_lateral_latency_ANOVA = [avg_lateral_latency_ANOVA; reshape(avg_lateral_latency,[],1)];
    avg_LR_latency_ANOVA = [avg_LR_latency_ANOVA; reshape(avg_LR_latency,[],1)];

    taste_lbls_ANOVA = [taste_lbls_ANOVA; reshape(repmat(tastes2use,6,1),[],1)];
    LR_lbls_ANOVA = [LR_lbls_ANOVA; reshape(repmat(["L","R"],6,1),[],1)];

    LR_for_ttest = [LR_for_ttest; avg_LR_latency];
end
anova1(avg_sampling_duration_ANOVA, taste_lbls_ANOVA)
%
anova1(avg_lateral_latency_ANOVA, taste_lbls_ANOVA)
%
% anova1(avg_LR_latency_ANOVA, LR_lbls_ANOVA)
[~, pval_latency] = ttest(LR_for_ttest(:,1), LR_for_ttest(:,2)) ;

% overall_lateral_latency
%% number of trials for each taste/animal
trial_counts_summary=[]; trial_counts_all=[]; trial_counts_avg=[];
fields_all = fieldnames(total_perf_all);
fields2use = string(fields_all(6:end));
for i=1:length(Animals)
    for j=1:length(fields2use)
        taste_labels = string(char(total_perf_all(i).(fields2use(j)).TasteID));
        [trial_counts, trial_lbls,~] = groupcounts(taste_labels);
        trial_lbls = strtrim(trial_lbls);
        for p=1:length(trial_lbls)
            trial_counts_summary(i).(trial_lbls(p)) = trial_counts(p);
        end
        trial_counts_all(i,:,j) = trial_counts;
        trial_counts_avg(i,j) = mean(trial_counts);
    end
end

%mean trial numbers reported in manuscript, pre and post learning
avg_prelearning = mean(trial_counts_avg(:,1:2),'all');
avg_prelearning = mean(trial_counts_avg(:,3:4),'all');

%range of trial numbers reported in manuscript, pre and post learning
range_prelearning = [min(trial_counts_all(:,:,1:2),[],'all') max(trial_counts_all(:,:,1:2),[],'all')];
range_postlearning = [min(trial_counts_all(:,:,3:4),[],'all') max(trial_counts_all(:,:,3:4),[],'all')];

