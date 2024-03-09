function plt_example_neuron_PSTH(neuron2use, time_stamps, neural_data, group2use, xlim2use)

    groups = string(fieldnames(neural_data));
    time_stamps_plt = mean(time_stamps.(groups(group2use)){1,neuron2use});    
    %start here
    data_plt1 = neural_data.(groups(group2use)){1,neuron2use};
    data_plt2 = neural_data.(groups(group2use)){2,neuron2use};
    
    plot_sem(time_stamps_plt, data_plt1)
    plot_sem(time_stamps_plt, data_plt2)
    p1=plot(time_stamps_plt, mean(data_plt1),'linewidth',2, 'Color','c');
    p2=plot(time_stamps_plt, mean(data_plt2),'linewidth',2, 'Color','k');

    box off
    axis square
    ylim2use=get(gca,'ylim');
    line([0 0],ylim2use, 'color','k','linestyle','--')

    legend([p1,p2],'Suc','NaCl','location','northwest')
    xticks2use = round(xlim2use);
    xticks(xticks2use(1):xticks2use(2));
    xlim(xlim2use)
    ylim(ylim2use)
%     yticks(-.1:.1:.3)
    ylabel('\DeltaF/f')
    xlabel('Time (sec)')
    title(['Example Neuron'])
end