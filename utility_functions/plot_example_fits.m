function plot_example_fits(neural_data, data2fit, f_stat, fit_vals, time_stamps,  neuron2use,  group2use, window2use_plt, resp_window, bin2use, subplts2use)
    groups = string(fieldnames(neural_data));
 
    line2 = {@(b,x) (b(1)*x./100) + b(2),... %linear function
                @(b,x) b(1)+(b(2)-b(1))./(1+10.^((50-x)*b(3))),...%sig function
                @(b,x) (b(1)*abs(b(1)*x-b(1)*50)/100)+b(2)}; %V shape
    lbls = ["Suc","75/25","60/40","40/60","25/75","0/100","Taste"];

    % neuron2use = 13;
    colors = cool(6);
    X = [100; 75; 60; 40; 25; 0]; %x values
    X2 = 100:-1:0;
    
    subplot(2,4,subplts2use(1)) %changed to 1 row
    hold on
    for j=1:6
        avg_activity=[];
        avg_activity = mean(neural_data.(groups(group2use)){j,neuron2use});
        ts_plt = mean(time_stamps.(groups(group2use)){j,neuron2use});

        plot(ts_plt, avg_activity,'color',colors(j,:),'linewidth',2)
    end
    xlim(window2use_plt)
    xticks(round(window2use_plt(1)):round(window2use_plt(2)));
    ylim2use = get(gca, 'ylim');
    line([0 0],ylim2use,'color','k','linestyle','--')
%     rectangle('Position',[resp_window(1) ylim2use(1) abs(diff(resp_window)) diff(ylim2use)],'facecolor',[.15 .15 .15 .1],'edgecolor','none')
    ylim(ylim2use)
    xlabel('Time (sec)')
    ylabel('\DeltaF/f')
    title('Example Neuron Activity')
    legend(lbls,'location','northeastoutside')
    axis square
    subplot(2,4,subplts2use(2))
    scatter(X, data2fit.(groups(group2use)){bin2use}(:,neuron2use),40,'filled','k')
    hold on
    plot(X2, line2{1}(fit_vals.(groups(group2use)){bin2use}{1}(:,neuron2use), X2),'color','k','linewidth',1)
    plot(X2, line2{2}(fit_vals.(groups(group2use)){bin2use}{2}(:,neuron2use), X2),'color','k','linewidth',1, 'linestyle','--')
    xticks(0:50:100)
    xlabel('Suc (mM)')
    ylabel('dff')
    title('Example Neuron Fits')
    axis square
    fstats2use = round(f_stat.(groups(group2use)){bin2use}(1:2,neuron2use),1);
    legend('Activity',['Linear fit, f = ', num2str(fstats2use(1))],...
        ['Sig fit, f = ', num2str(fstats2use(2))],'location','northeastoutside')
end