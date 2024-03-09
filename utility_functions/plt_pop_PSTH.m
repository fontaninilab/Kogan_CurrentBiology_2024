function plt_pop_PSTH(neurons2use, time_stamps, neural_data, group2use, xlim2use,delta)
    groups = string(fieldnames(neural_data));
    for i=1:length(neurons2use)
       pop_PSTH(i,:) = mean(neural_data.(groups(group2use)){1,neurons2use(i)});
       pop_PSTH_all_trials(i,:) = mean([neural_data.(groups(group2use)){1,neurons2use(i)}; neural_data.(groups(group2use)){2,neurons2use(i)}]);
       delta_PSTH(i,:) = abs(mean(neural_data.(groups(group2use)){1,neurons2use(i)}) - mean(neural_data.(groups(group2use)){2,neurons2use(i)}));
       pop_PSTH_ts(i,:) = mean(time_stamps.(groups(group2use)){1,neurons2use(i)});
    end
    if delta==1
        plot(mean(pop_PSTH_ts), mean(delta_PSTH),'linewidth',1','color','k')
        plot_sem(mean(pop_PSTH_ts), delta_PSTH)
    else
        plot(mean(pop_PSTH_ts), mean(pop_PSTH_all_trials),'linewidth',1','color','k')
        plot_sem(mean(pop_PSTH_ts), pop_PSTH_all_trials)
    end
    % rectangle('Position',[0 ylim2use(1) 1.5 diff(ylim2use)],'facecolor',[.15 .15 .15 .1],'edgecolor','none')
    ylim2use = get(gca, 'ylim');
    line([0 0], ylim2use, 'linestyle','--','color','k')
    box off
    axis square
    xlabel('Time (sec)')
    ylabel('\DeltaF/f')
    xlim(xlim2use)
    xticks2use = round(xlim2use);
    xticks(xticks2use(1):xticks2use(2))
    ylim(ylim2use)
%     yticks(0:0.05:0.15)
%     title(['Sampling Responsive',newline,'Pop PSTH'])
    title('Pop PSTH')
end